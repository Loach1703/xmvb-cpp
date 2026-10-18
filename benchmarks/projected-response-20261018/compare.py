#!/usr/bin/env python3
"""Compare native projected-response runs with saved 66a31968 results."""

import argparse
import csv
from pathlib import Path


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path, nargs="+")
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    refs = {(row["case"], row["mode"]): row
            for name in ("ri.csv", "exact.csv")
            for row in read(args.baseline / name)}
    rows = []
    for path in args.candidate:
        for new in read(path):
            old = refs[new["case"], new["mode"]]
            if old["input_sha256"] != new["input_sha256"]:
                raise RuntimeError(f"input mismatch: {new['case']}")
            row = dict(case=new["case"], mode=new["mode"],
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
            print(f"{new['case']} {new['mode']}: {old['process_s']} -> "
                  f"{new['process_s']} s; {row['speedup']:.3f}x; "
                  f"{old['steps']}/{old['hvp']} -> {new['steps']}/{new['hvp']}; "
                  f"dE={row['delta_Eh']:.3g} Eh; matched={row['same_node_mask']}")
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
