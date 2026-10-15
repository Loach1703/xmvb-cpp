#!/usr/bin/env python3
"""Compare native CIAH cleanup runs with saved runs, without rerunning baselines."""
import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CASES = {"241": "241_vb", "MnF2": "MnF2", "FeCl2": "FeCl2",
         "QOCGUM": "csd_QOCGUM_vb", "CERRAS": "csd_CERRAS_vb",
         "LOFLEA": "csd_LOFLEA_vb"}


def read(path):
    return path.read_text() if path.exists() else ""


def capture(pattern, text):
    match = re.search(pattern, text)
    return match.group(1) if match else ""


def native(folder):
    text = read(folder / "output.out")
    trace = list(csv.DictReader(read(folder / "neo.tsv").splitlines(), delimiter="\t"))
    reason = capture(r"Termination reason\s*:\s*([^\n]+)", text)
    status = "running" if folder.exists() else "not_submitted"
    if (folder / "run.exit").exists():
        status = "converged" if read(folder / "run.exit").strip() == "0" and reason in {
            "neo_projected_gradient_tolerance", "neo_initial_tolerance"} else "failed"
    rss = capture(r"Maximum resident set size \(kbytes\):\s*(\d+)", read(folder / "resource.time"))
    return dict(status=status, steps=len(trace),
                hvp=sum(int(r["orbital_hvp_actions"]) for r in trace),
                energy=capture(r"Final total energy\s*:\s*(-?[\d.]+)", text),
                wall_s=capture(r"End-to-end wall time\s*:\s*([\d.]+)", text),
                rss_MiB=float(rss)/1024 if rss else "")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("current")
    parser.add_argument("--previous", default="riesz-reuse-1250148")
    parser.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES))
    args = parser.parse_args()
    current = ROOT / args.current
    history = {r["case"]: r for r in csv.DictReader(
        (ROOT.parents[1] / "history/paired_runs.csv").open())
        if r["method"] == "neo" and r["revision"] == "c2776c4"}
    rows = []
    for case in args.cases:
        alias = CASES[case]
        row = {"case": case, **native(current / case)}
        previous = native(ROOT / args.previous / case)
        row.update({"previous_" + key: value for key, value in previous.items()})
        row["energy_delta_Eh"] = (float(row["energy"]) - float(previous["energy"])) \
            if row["energy"] and previous["energy"] else ""
        old = history.get(alias, {})
        row["historical_hvp"] = old.get("orbital_hvp_actions", "")
        row["historical_wall_s"] = old.get("total_wall_s", "")
        row["historical_energy"] = old.get("energy_final_hartree", "")
        row["historical_steps"] = old.get("accepted_steps", "")
        row["historical_rss_MiB"] = float(old["peak_rss_kb"])/1024 \
            if old.get("peak_rss_kb") else ""
        row["note"] = "Historical Hanhai25; native runs Hanhai22, not matched nodes. FeCl2 historical endpoint differs."
        rows.append(row)
    writer = csv.DictWriter(sys.stdout, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
