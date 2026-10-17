#!/usr/bin/env python3
"""Compare saved, matched GitHub/current CAS12 native runs; never run binaries."""

import argparse
import csv
from pathlib import Path
import re

import numpy as np


BLOCKS = (
    "action_hamiltonian", "action_overlap", "accepted_gradient_overlap",
    "accepted_gradient_one_electron", "accepted_gradient_pair_kernel",
    "direction_hamiltonian", "direction_overlap", "reverse_overlap",
    "reverse_one_electron", "reverse_pair_kernel",
)


def fields(path):
    return dict(line.split("=", 1) for line in path.read_text().splitlines() if "=" in line)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    commit = (args.run / "baseline-commit.txt").read_text().strip()
    rows = []
    for count in (1, 3, 20, 132):
        info, data, rss = {}, {}, {}
        for variant in ("baseline", "current"):
            root = args.run / variant
            info[variant] = fields(root / f"{count}.out")
            data[variant] = np.fromfile(root / f"{count}.bin", dtype="<f8")
            resource = (root / f"{count}.time").read_text()
            rss[variant] = int(re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", resource)[1]) / 1024
        old, new = info["baseline"], info["current"]
        if len(data["baseline"]) != len(data["current"]):
            raise ValueError(f"different dump sizes at {count} structures")
        row = dict(structures=count, baseline_commit=commit,
                   node=(args.run / "node.txt").read_text().strip(), actual_threads=8,
                   scope="constructed_CAS12_fixed_point_not_molecular_SCF",
                   baseline_s=float(old["total_seconds"]), current_s=float(new["total_seconds"]),
                   baseline_rss_mib=rss["baseline"], current_rss_mib=rss["current"])
        row["wall_speedup"] = row["baseline_s"] / row["current_s"]
        for stage in ("construction", "action", "accepted_prepare_gradient",
                      "orbital_direction", "cached_reverse"):
            for variant in ("baseline", "current"):
                row[f"{variant}_{stage}_s"] = float(info[variant][stage + "_seconds"])
        for variant in ("baseline", "current"):
            row[f"{variant}_direction_reverse_s"] = (
                row[f"{variant}_orbital_direction_s"] + row[f"{variant}_cached_reverse_s"])
        offset = 0
        largest_relative = 0.0
        for block in BLOCKS:
            size = int(old[f"array.{block}.count"])
            for key in ("count", "logical_rows", "logical_columns", "layout"):
                if old[f"array.{block}.{key}"] != new[f"array.{block}.{key}"]:
                    raise ValueError(f"different {block} shape/layout at {count}")
            a, b = (data[v][offset:offset + size] for v in ("baseline", "current"))
            if not (np.isfinite(a).all() and np.isfinite(b).all()):
                raise ValueError(f"nonfinite {block} at {count}")
            row[f"{block}_max_abs_error"] = float(np.max(np.abs(a - b), initial=0.0))
            largest_relative = max(largest_relative, float(np.linalg.norm(a - b)) /
                                   max(1.0, float(np.linalg.norm(a)), float(np.linalg.norm(b))))
            offset += size
        if offset != len(data["baseline"]):
            raise ValueError("unaccounted numerical dump payload")
        row["largest_block_relative_error"] = largest_relative
        row["parity_passed"] = largest_relative <= 3e-10
        rows.append(row)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(row["structures"], "total", row["baseline_s"], row["current_s"],
              "RSS", row["baseline_rss_mib"], row["current_rss_mib"],
              "relative_error", row["largest_block_relative_error"])
    if not all(row["parity_passed"] for row in rows):
        raise SystemExit("GitHub/current numerical parity failed; retain failed rows")


if __name__ == "__main__":
    main()
