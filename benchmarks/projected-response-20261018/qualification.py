#!/usr/bin/env python3
"""Summarize saved Slurm qualification outputs without executing C++ code."""

import argparse
import csv
import re
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    ctest = (args.artifact / "ctest.out").read_text()
    summary = re.search(r"100% tests passed, 0 tests failed out of (\d+)", ctest)
    if not summary:
        raise RuntimeError("qualification test suite did not pass")
    source_check = (args.artifact / "source-check.txt").read_text().splitlines()
    manifest = (args.artifact / "source.sha256").read_text().splitlines()
    if len(source_check) != len(manifest) or not all(line.endswith(": OK") for line in source_check):
        raise RuntimeError("source manifest verification is incomplete")
    keys = ("gradient_energy_fd_absolute_error", "max_rel_diff", "core_max_rel_diff",
            "outer_max_rel_diff", "component_decomposition_error",
            "coupling_batch_orbital_error", "coupling_batch_structure_error",
            "orbital_coupling_symmetry_error", "relaxed_hvp_symmetry_error")
    rows = []
    for mode in ("exact", "ri"):
        for case in ("F2", "FeCl2"):
            values = dict(line.split(" = ", 1) for line in
                          (args.artifact / f"{mode}-{case}" / "hvp.out").read_text().splitlines()
                          if " = " in line)
            row = dict(mode=mode, case=case, ctest_passed=int(summary[1]),
                       source_files_verified=len(source_check))
            row.update({key: values[key] for key in keys})
            if float(row["max_rel_diff"]) > 1e-5:
                raise RuntimeError(f"HVP finite difference failed: {case}/{mode}")
            rows.append(row)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
