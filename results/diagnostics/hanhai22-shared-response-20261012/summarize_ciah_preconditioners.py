#!/usr/bin/env python3
"""Print paired native CIAH run metrics as CSV; missing data stay empty."""

import argparse
import csv
from pathlib import Path
import re
import sys


def read(path):
    return path.read_text(errors="replace") if path.exists() else ""


def capture(pattern, text):
    match = re.search(pattern, text)
    return match[1] if match else ""


def seconds(value):
    if not value:
        return ""
    total = 0.0
    for part in value.split(":"):
        total = 60 * total + float(part)
    return total


def row(system, preconditioner, path):
    output = read(path / "output.out")
    resources = read(path / "resource.time")
    exit_code = read(path / "run.exit").strip()
    trace = path / "neo.tsv"
    records = list(csv.DictReader(trace.open(), delimiter="\t")) if trace.exists() else []
    rss = capture(r"Maximum resident set size \(kbytes\):\s*(\d+)", resources)
    return {
        "system": system,
        "preconditioner": preconditioner,
        "status": "finished" if exit_code else "running" if path.exists() else "pending",
        "converged": bool(re.search(r"VBSCF converged in\s+\d+ iterations", output)),
        "steps": capture(r"VBSCF converged in\s+(\d+) iterations", output),
        "recorded_steps": len(records) if records else "",
        "recorded_orbital_hvp": sum(int(record["orbital_hvp_actions"]) for record in records) if records else "",
        "energy_Eh": capture(r"Final total energy\s*:\s*(-?[\d.]+)", output),
        "wall_s": seconds(capture(r"Elapsed \(wall clock\).*\):\s*(\S+)", resources)),
        "peak_RSS_MiB": float(rss) / 1024 if rss else "",
        "exit_code": exit_code,
        "termination": capture(r"Termination reason\s*:\s*([^\n]+)", output).strip(),
        "node": (read(path / "node.txt") or read(path.parent / "node.txt")).strip(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    args = parser.parse_args()
    rows = [row(system, preconditioner, args.run_root / system / preconditioner)
            for system in ("F2", "241", "240", "MnF2", "FeCl2")
            for preconditioner in ("hessian-diagonal", "one-electron")]
    writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
