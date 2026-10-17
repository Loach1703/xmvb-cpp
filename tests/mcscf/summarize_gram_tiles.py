#!/usr/bin/env python3
"""Parse saved source-factor kernel benchmarks; does not launch workloads."""

import argparse
import csv
import sys
from pathlib import Path


def summarize(root):
    node = (root / "node.txt").read_text().strip()
    for count in (1, 3, 20, 132):
        fields = dict(
            line.split("=", 1)
            for line in (root / f"{count}.out").read_text().splitlines()
            if "=" in line
        )
        yield {
            "run": root.name,
            "node": node,
            "structures": int(fields["structures"]),
            "touched_products": int(fields["products"]),
            **{
                key: float(fields[key])
                for key in (
                    "gram_action_seconds", "gram_tangent_seconds",
                    "gram_adjoint_seconds", "gram_curvature_seconds",
                    "gram_duality_error", "gram_curvature_fd_error",
                )
            },
            "rss_mib": int(fields["max_rss_kib"]) / 1024,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path, nargs="+")
    args = parser.parse_args()
    rows = [row for root in args.runs for row in summarize(root)]
    writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
