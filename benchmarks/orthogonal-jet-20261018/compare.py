#!/usr/bin/env python3
"""Compare candidate native runs against existing, input-matched records."""

import argparse
import csv
from pathlib import Path


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", nargs="+", type=Path)
    parser.add_argument("--baseline", action="append", type=Path, required=True,
                        help="Saved CSVs in priority order, most recent last")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    refs = {}
    for path in args.baseline:
        for old in read(path):
            refs[old["case"], old["mode"]] = old, str(path)
    rows = []
    for path in args.candidate:
        for new in read(path):
            old, source = refs[new["case"], new["mode"]]
            if new["input_sha256"] != old["input_sha256"]:
                raise RuntimeError(f"input mismatch: {new['case']}")
            row = dict(case=new["case"], mode=new["mode"], baseline_csv=source,
                       baseline_job=old["job"], candidate_job=new["job"],
                       baseline_node=old["node"], candidate_node=new["node"],
                       same_node_mask=old["node"] == new["node"] and
                           old["affinity"] == new["affinity"])
            for key in ("steps", "hvp", "process_s", "rss_GiB", "final_Eh",
                        "ci_actions", "ci_s", "hvp_s", "effective_cores"):
                row["baseline_" + key] = old[key]
                row["candidate_" + key] = new[key]
            row["speedup"] = float(old["process_s"]) / float(new["process_s"])
            row["delta_Eh"] = float(new["final_Eh"]) - float(old["final_Eh"])
            row["converged"] = new["converged"]
            rows.append(row)
            print(f"{new['case']}/{new['mode']}: {old['process_s']} -> "
                  f"{new['process_s']} s; {row['speedup']:.3f}x; "
                  f"{old['steps']}/{old['hvp']} -> {new['steps']}/{new['hvp']}; "
                  f"dE={row['delta_Eh']:.3g} Eh")
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
