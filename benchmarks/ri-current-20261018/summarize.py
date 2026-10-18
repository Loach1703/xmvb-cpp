#!/usr/bin/env python3
"""Read native RI logs and saved exact references; never launch a calculation."""

import argparse
import csv
from pathlib import Path


def fields(path):
    return dict(line.strip().split(": ", 1) for line in path.read_text().splitlines()
                if ": " in line)


def elapsed(value):
    seconds = 0.0
    for field in value.split(":"):
        seconds = seconds * 60 + float(field)
    return seconds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path)
    parser.add_argument("references", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with args.references.open() as stream:
        references = {r["case"]: r for r in csv.DictReader(stream)}
    rows = []
    for case in ("241", "c12-cov", "c14-cov", "FeCl2", "MnF2"):
        path = args.artifact / case
        report = {k.strip(): v.strip() for k, v in fields(path / "report-tail.txt").items()}
        resource = fields(path / "resource.time")
        hashes = {name: sha for sha, name in
                  (line.split() for line in (path / "hashes.txt").read_text().splitlines())}
        with (path / "ciah.tsv").open() as stream:
            trace = list(csv.DictReader(stream, delimiter="\t"))
        process_s = elapsed(resource["Elapsed (wall clock) time (h:mm:ss or m:ss)"])
        user_s, system_s = (float(resource[key]) for key in
                            ("User time (seconds)", "System time (seconds)"))
        converged = (path / "run.exit").read_text().strip() == "0" and \
            report.get("Termination reason") == "ciah_projected_gradient_tolerance"
        counters = dict(ci_actions=0, ci_s=0.0, hvp_s=0.0, keyframe_s=0.0)
        for line in (path / "run.err").read_text().splitlines():
            columns = line.split(",")
            if len(columns) < 5 or columns[0] != "CIAH_DIAG":
                continue
            values = dict(item.split("=", 1) for item in columns[4:] if "=" in item)
            if columns[3] == "ciah_ci_solved":
                counters["ci_actions"] += int(values["actions"])
                counters["ci_s"] += float(values["seconds"])
            elif columns[3] == "ciah_expansion_cost":
                counters["hvp_s"] += float(values["hvp_seconds"])
            elif columns[3] == "ciah_keyframe_cost":
                counters["keyframe_s"] += float(values["seconds"])
        source = (path / "input.xmi").read_text().splitlines()
        basis = next(line.split("=", 1)[1].strip() for line in source
                     if line.strip().lower().startswith("basis="))
        row = dict(case=case, mode="ri", job=(path / "slurm-job.txt").read_text().strip(),
                   node=(path / "node.txt").read_text().strip(),
                   affinity=(path / "affinity.txt").read_text().split(":", 1)[1].strip(),
                   basis=basis, auxbasis=basis.lower()+"-jkfit", converged=converged,
                   steps=len(trace), hvp=sum(int(t["orbital_hvp_actions"]) for t in trace),
                   per_step_hvp=":".join(t["orbital_hvp_actions"] for t in trace),
                   initial_Eh=report.get("Initial total energy", ""),
                   final_Eh=report.get("Final total energy", ""),
                   process_s=process_s, application_s=report["End-to-end wall time"].split()[0],
                   scf_s=report["SCF iteration wall time"].split()[0],
                   rss_GiB=int(resource["Maximum resident set size (kbytes)"])/1048576,
                   effective_cores=(user_s+system_s)/process_s, user_s=user_s, system_s=system_s,
                   input_sha256=hashes["input.xmi"],
                   binary_sha256=next(sha for name, sha in hashes.items()
                                      if name.endswith("/xmvb-cpp.exe")), **counters)
        reference = references[case]
        if row["input_sha256"] != reference["input_sha256"]:
            raise RuntimeError(f"{case}: input differs from saved exact reference")
        for key in ("job", "node", "source", "process_s", "rss_GiB", "steps", "hvp", "final_Eh"):
            row[f"exact_{key}"] = reference[key]
        row["reference_over_ri_time"] = float(reference["process_s"])/process_s if converged else ""
        row["delta_E_mEh"] = 1000*(float(row["final_Eh"])-float(reference["final_Eh"])) if converged else ""
        rows.append(row)
        print(f"{case}: converged={converged}, {row['steps']}/{row['hvp']}, "
              f"{process_s:.2f}s, {row['rss_GiB']:.3f}GiB, "
              f"mean cores={row['effective_cores']:.2f}, delta={row['delta_E_mEh']}mEh")
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
