#!/usr/bin/env python3
"""Summarize saved selected-space H/S kernels, keeping native and parity runs apart."""

import argparse
import csv
import sys
from pathlib import Path


def fields(path):
    return dict(line.split("=", 1) for line in path.read_text().splitlines() if "=" in line)


def summarize(root):
    for count in (1, 3, 20, 132):
        native_path = root / f"{count}.out"
        reference_path = root / f"{count}-parity.out"
        if not native_path.exists():
            native_path = root / f"{count}-selected12.out"
            reference_path = root / f"{count}-selected12-parity.out"
        native = fields(native_path)
        reference = fields(reference_path)
        yield {
            "run": root.name,
            "node": (root / "node.txt").read_text().strip(),
            "structures": count,
            "cold_seconds": float(native["cold_apply_seconds"]),
            "warm_seconds": float(native["warm_apply_seconds"]),
            "rss_mib": int(native["max_rss_kib"]) / 1024,
            "rayleigh_energy": float(native["rayleigh_energy"]),
            "parity_process_warm_seconds": float(reference["warm_apply_seconds"]),
            "hamiltonian_error": float(reference["hamiltonian_reference_error"]),
            "overlap_error": float(reference["overlap_reference_error"]),
            "source_columns": int(native.get("source_columns", 0)),
            "source_tiles": int(native.get("source_tiles", 0)),
            "source_workspace_bound": int(native.get("source_workspace_bound", 0)),
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
