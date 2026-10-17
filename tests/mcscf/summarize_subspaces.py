#!/usr/bin/env python3
"""Summarize saved 240 subspace runs; never infer convergence from job status."""

import argparse
import csv
import re
import sys
from pathlib import Path


def field(text, pattern):
    match = re.search(pattern, text)
    if match is None:
        raise ValueError(f"Missing result field: {pattern}")
    return match.group(1)


def summarize(root):
    for label in ("cov", "ion01", "ion02", "ion03", "full"):
        run = root / label
        output = (run / "output.out").read_text()
        timing = (run / "resource.time").read_text()
        with (run / "neo.tsv").open() as stream:
            trace = list(csv.DictReader(stream, delimiter="\t"))
        elapsed = field(timing, r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)")
        seconds = 0.0
        for part in elapsed.split(":"):
            seconds = seconds * 60 + float(part)
        yield {
            "space": label,
            "structures": int(field(output, r"Number of structures:\s+(\d+)")),
            "converged": bool(re.search(r"VBSCF converged in\s+\d+", output)),
            "iterations": field(output, r"VBSCF converged in\s+(\d+)"),
            "orbital_hvp": sum(int(row["orbital_hvp_actions"]) for row in trace),
            "initial_energy_Eh": field(output, r"Initial total energy\s*:\s*([-\d.]+)"),
            "final_energy_Eh": field(output, r"Final total energy\s*:\s*([-\d.]+)"),
            "scf_seconds": field(output, r"SCF iteration wall time\s*:\s*([\d.]+)"),
            "application_seconds": field(output, r"End-to-end wall time\s*:\s*([\d.]+)"),
            "process_seconds": seconds,
            "rss_mib": int(field(timing, r"Maximum resident set size \(kbytes\): (\d+)")) / 1024,
            "node": (run / "node.txt").read_text().strip(),
            "exit_code": (run / "run.exit").read_text().strip(),
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    rows = list(summarize(args.directory))
    writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
