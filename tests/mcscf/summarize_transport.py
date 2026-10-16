#!/usr/bin/env python3
"""Aggregate saved scalar CI-transport records, without launching calculations."""

import argparse
import csv
from pathlib import Path
import sys


def summarize(root):
    node = (root / "node.txt").read_text().strip()
    for directory in (root / "native", root / "phases"):
        if not directory.is_dir():
            continue
        for case in sorted(directory.iterdir()):
            log = case / "run.err"
            if not log.is_file():
                continue
            records = [dict(field.split("=", 1) for field in line.split(",")[1:])
                       for line in log.read_text().splitlines()
                       if line.startswith("XMVB_TRANSPORT_SOLVE,")]
            if not records:
                continue
            yield dict(
                run=root.name, mode=directory.name, space=case.name, node=node,
                calls=len(records),
                metric_actions=sum(int(row["actions"]) for row in records),
                cg_iterations=sum(int(row["iterations"]) for row in records),
                certifications=sum(int(row["certifications"]) for row in records),
                max_relative_residual=max(float(row["relative_residual"]) for row in records),
                solve_seconds=sum(float(row["wall"]) for row in records),
                metric_action_seconds=sum(float(row["action_wall"]) for row in records))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path, nargs="+")
    args = parser.parse_args()
    rows = [row for root in args.runs for row in summarize(root)]
    if rows:
        writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
