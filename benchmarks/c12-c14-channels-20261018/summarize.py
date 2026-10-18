#!/usr/bin/env python3
"""Join native C12/C14 results to frozen records without rerunning a baseline."""

import argparse
import csv
from collections import deque
from pathlib import Path


def elapsed(value):
    result = 0.0
    for field in value.split(":"):
        result = 60.0 * result + float(field)
    return result


def fields(path, tail=None):
    with path.open() as stream:
        lines = deque(stream, maxlen=tail) if tail else stream
        return dict(line.strip().split(": ", 1) for line in lines if ": " in line)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("job")
    parser.add_argument("baseline", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with args.baseline.open() as stream:
        baseline = list(csv.DictReader(stream))
    current, comparisons = [], []
    for old in baseline:
        case = old["case"]
        path = args.root / f"run-{args.job}" / case
        if (path / "run.exit").read_text().strip() != "0":
            raise RuntimeError(f"{case}: application failed")
        report_path = path / "report-tail.txt"
        if not report_path.exists():
            report_path = path / "output.out"
        report = fields(report_path, tail=90)
        report = {k.strip(): v.strip() for k, v in report.items()}
        if report["Termination reason"] != "ciah_projected_gradient_tolerance":
            raise RuntimeError(f"{case}: not converged")
        resource = fields(path / "resource.time")
        hashes = {name: sha for sha, name in
                  (line.split() for line in (path / "hashes.txt").read_text().splitlines())}
        if hashes["input.xmi"] != old["input_sha256"]:
            raise RuntimeError(f"{case}: input differs from the baseline")
        with (path / "ciah.tsv").open() as stream:
            trace = list(csv.DictReader(stream, delimiter="\t"))
        ci_actions, ci_seconds, hvp_seconds, keyframe_seconds = 0, 0.0, 0.0, 0.0
        for line in (path / "run.err").read_text().splitlines():
            columns = line.split(",")
            if len(columns) < 5 or columns[0] != "CIAH_DIAG":
                continue
            values = dict(item.split("=", 1) for item in columns[4:] if "=" in item)
            if columns[3] == "ciah_ci_solved":
                ci_actions += int(values["actions"])
                ci_seconds += float(values["seconds"])
            elif columns[3] == "ciah_expansion_cost":
                hvp_seconds += float(values["hvp_seconds"])
            elif columns[3] == "ciah_keyframe_cost":
                keyframe_seconds += float(values["seconds"])
        row = dict(case=case, job=f"{args.job}_{len(current)}",
                   node=(path / "node.txt").read_text().strip(),
                   steps=len(trace), hvp=sum(int(t["orbital_hvp_actions"]) for t in trace),
                   per_step_hvp=":".join(t["orbital_hvp_actions"] for t in trace),
                   initial_Eh=float(report["Initial total energy"]),
                   final_Eh=float(report["Final total energy"]),
                   application_s=float(report["End-to-end wall time"].split()[0]),
                   scf_s=float(report["SCF iteration wall time"].split()[0]),
                   process_s=elapsed(resource["Elapsed (wall clock) time (h:mm:ss or m:ss)"]),
                   rss_GiB=int(resource["Maximum resident set size (kbytes)"]) / 1048576.0,
                   user_s=float(resource["User time (seconds)"]),
                   system_s=float(resource["System time (seconds)"]),
                   ci_actions=ci_actions, ci_s=ci_seconds, hvp_s=hvp_seconds,
                   keyframe_s=keyframe_seconds, input_sha256=hashes["input.xmi"],
                   binary_sha256=hashes["bin/xmvb-cpp.exe"])
        if abs(row["initial_Eh"] - float(old["initial_Eh"])) > 1e-10 or \
                abs(row["final_Eh"] - float(old["final_Eh"])) > 5e-10:
            raise RuntimeError(f"{case}: energy regression")
        current.append(row)
        comparisons.append(dict(case=case, before_job=old["job"], after_job=row["job"],
            **{f"{prefix}_{key}": source[key]
               for prefix, source in (("before", old), ("after", row))
               for key in ("node", "steps", "hvp", "process_s", "application_s", "rss_GiB", "final_Eh")},
            speedup=float(old["process_s"]) / row["process_s"],
            delta_Eh=row["final_Eh"] - float(old["final_Eh"])))
    for path, rows in ((args.output, current),
                       (args.output.with_name(args.output.stem + "-comparison.csv"), comparisons)):
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    for row in comparisons:
        print(row)


if __name__ == "__main__":
    main()
