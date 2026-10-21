#!/usr/bin/env python3
"""Compare saved native RI runs without launching or repeating workloads."""

import argparse
import csv
from pathlib import Path


def read(path):
    with path.open() as stream:
        return {row["case"]: row for row in csv.DictReader(stream)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    baseline, candidate = read(args.baseline), read(args.candidate)
    rows = []
    for case, current in candidate.items():
        old = baseline[case]
        for key in ("input_sha256", "basis", "auxbasis", "mode"):
            if old[key] != current[key]:
                raise RuntimeError(f"{case}: different {key}")
        if old["converged"] != "True" or current["converged"] != "True":
            raise RuntimeError(f"{case}: comparison requires converged runs")
        row = {"case": case}
        for name, run in (("baseline", old), ("candidate", current)):
            for key in ("job", "node", "affinity", "binary_sha256", "steps", "hvp",
                        "process_s", "rss_GiB", "effective_cores", "final_Eh",
                        "ci_actions", "ci_s", "hvp_s", "keyframe_s"):
                row[f"{name}_{key}"] = run[key]
        row["same_node"] = old["node"] == current["node"]
        row["wall_speedup"] = float(old["process_s"])/float(current["process_s"])
        row["rss_ratio"] = float(current["rss_GiB"])/float(old["rss_GiB"])
        row["energy_difference_Eh"] = float(current["final_Eh"])-float(old["final_Eh"])
        rows.append(row)
        print(f"{case}: {old['steps']}/{old['hvp']} -> "
              f"{current['steps']}/{current['hvp']}, "
              f"{old['process_s']} -> {current['process_s']} s, "
              f"speedup={row['wall_speedup']:.3f}, same_node={row['same_node']}")
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
