#!/usr/bin/env python3
"""Read frozen CIAH tail experiments; write CSV to stdout, never rerun a case."""

import csv
import io
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[2]
RUNS = (
    ("retained", "291714", "ciah-path-20261017/borrow"),
    ("one-electron", "291729", "ciah-tail-20261017/one-electron"),
    ("current-seed", "291751", "ciah-tail-20261017/current-seed"),
    ("base-only", "291758", "ciah-tail-20261017/base-only"),
)
CASES = ("240", "241", "FeCl2", "F2", "MnF2")
FIELDS = (
    "case", "variant", "job", "node", "steps", "hvp", "keyframes",
    "ritz_solves", "final_energy_Eh", "wall_seconds", "peak_rss_KiB",
    "input_sha256", "binary_sha256", "artifact_path",
)


def value(pattern, text):
    match = re.search(pattern, text, flags=re.MULTILINE)
    if match is None:
        raise ValueError(f"Missing field: {pattern}")
    return match[1]


def summarize(directory):
    output = (directory / "output.out").read_text()
    if (directory / "run.exit").read_text().strip() != "0":
        raise ValueError(f"Failed run: {directory}")
    value(r"^Termination reason\s*:\s*(ciah_projected_gradient_tolerance)$", output)
    rows = list(csv.DictReader(io.StringIO(
        (directory / "ciah.tsv").read_text().replace("\\t", "\t")), delimiter="\t"))
    if not rows:
        raise ValueError(f"Missing trajectory: {directory}")
    resource = (directory / "resource.time").read_text()
    elapsed = value(r"Elapsed \(wall clock\) time.*:\s*(\S+)\s*$", resource)
    seconds = 0.0
    for part in elapsed.split(":"):
        seconds = 60.0 * seconds + float(part)
    hashes = dict(line.split(maxsplit=1)[::-1]
                  for line in (directory / "hashes.txt").read_text().splitlines())
    binary = next(digest for name, digest in hashes.items() if name.endswith("/xmvb-cpp.exe"))
    return {
        "node": (directory / "node.txt").read_text().strip(),
        "steps": len(rows),
        "hvp": sum(int(row["orbital_hvp_actions"]) for row in rows),
        "keyframes": sum(int(row["keyframes"]) for row in rows),
        "ritz_solves": sum(int(row["micro_iterations"]) for row in rows),
        "final_energy_Eh": value(r"^Final total energy\s*:\s*([-\d.]+)$", output),
        "wall_seconds": f"{seconds:.2f}",
        "peak_rss_KiB": value(r"Maximum resident set size \(kbytes\):\s*(\d+)", resource),
        "input_sha256": hashes["input.xmi"],
        "binary_sha256": binary,
        "artifact_path": str(directory.relative_to(ROOT)),
    }


def main():
    writer = csv.DictWriter(sys.stdout, fieldnames=FIELDS)
    writer.writeheader()
    inputs = {}
    for variant, job, path in RUNS:
        for index, case in enumerate(CASES):
            row = summarize(ROOT / "results/diagnostics" / path / case)
            if row["input_sha256"] != inputs.setdefault(case, row["input_sha256"]):
                raise ValueError(f"Unmatched input: {variant}/{case}")
            writer.writerow(dict(case=case, variant=variant, job=f"{job}_{index}", **row))


if __name__ == "__main__":
    main()
