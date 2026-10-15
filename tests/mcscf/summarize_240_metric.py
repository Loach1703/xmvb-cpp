#!/usr/bin/env python3
"""Summarize saved 240 runs; never launches calculations or old executables."""

import argparse
import csv
import re
import sys
from pathlib import Path


def match(text, pattern):
    found = re.search(pattern, text, re.MULTILINE)
    if not found:
        raise ValueError(f"Missing output field: {pattern}")
    return found.group(1)


def summarize(root):
    node = (root / "node.txt").read_text().strip()
    for space in ("cov", "ion01", "ion02", "ion03", "full"):
        case = root / "native" / space
        if not (case / "run.exit").exists():
            print(f"Not completed: {case}", file=sys.stderr)
            continue
        output = (case / "output.out").read_text()
        resource = (case / "resource.time").read_text()
        with (case / "neo.tsv").open() as stream:
            trace = list(csv.DictReader(stream, delimiter="\t"))
        elapsed = match(resource, r"Elapsed \(wall clock\).*?: ([0-9:.]+)")
        seconds = 0.0
        for field in elapsed.split(":"):
            seconds = 60 * seconds + float(field)
        yield {
            "run": root.name,
            "space": space,
            "node": node,
            "structures": int(match(output, r"Number of structures:\s+(\d+)")),
            "iterations": int(match(output, r"VBSCF converged in\s+(\d+) iterations")),
            "orbital_hvp": sum(int(row["orbital_hvp_actions"]) for row in trace),
            "initial_energy_Eh": float(match(output, r"Initial total energy\s+:\s+(-?[0-9.]+)")),
            "final_energy_Eh": float(match(output, r"Final total energy\s+:\s+(-?[0-9.]+)")),
            "scf_seconds": float(match(output, r"SCF iteration wall time\s+:\s+([0-9.]+)")),
            "application_seconds": float(match(output, r"End-to-end wall time\s+:\s+([0-9.]+)")),
            "process_seconds": seconds,
            "rss_mib": int(match(resource, r"Maximum resident set size \(kbytes\): (\d+)")) / 1024,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path, nargs="+")
    args = parser.parse_args()
    rows = [row for root in args.runs for row in summarize(root)]
    if rows:
        writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
