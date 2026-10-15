#!/usr/bin/env python3
"""Summarize native source-cover actions; reference-process RSS stays separate."""

import argparse
import csv
from pathlib import Path


def values(path):
    return dict(line.strip().split("=", 1) for line in path.read_text().splitlines() if "=" in line)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with args.baseline.open() as handle:
        baseline = {int(row["structures"]): row for row in csv.DictReader(handle)
                    if row["run"] == "source-tiles-1260459"}
    cases = [(f"selected12-{n}", f"{n}-selected12.out", n) for n in (1, 3, 20, 132)]
    cases += [("covalent14", "14-source-size.out", None)]
    cases += [(f"cross{n}", f"{n}-source-cover-size.out", None) for n in (12, 14)]
    rows = []
    for label, filename, count in cases:
        data = values(args.artifact / filename)
        row = dict(case=label, node=(args.artifact / "node.txt").read_text().strip(),
                   warm_ms=1000 * float(data["warm_apply_seconds"]),
                   rss_mib=float(data["max_rss_kib"]) / 1024,
                   source_factors=int(data["source_columns"]),
                   source_workspace_mib=float(data["source_workspace_bound"]) / 1024**2,
                   energy_hartree=float(data["rayleigh_energy"]),
                   previous_warm_ms="", previous_rss_mib="", energy_delta_hartree="",
                   h_error="", m_error="")
        if count is not None:
            old = baseline[count]
            row.update(previous_warm_ms=1000 * float(old["warm_seconds"]),
                       previous_rss_mib=float(old["rss_mib"]),
                       energy_delta_hartree=float(data["rayleigh_energy"]) - float(old["rayleigh_energy"]))
        if label.startswith("cross"):
            reference = values(args.artifact / filename.replace("-size", "-reference"))
            row.update(h_error=float(reference["hamiltonian_reference_error"]),
                       m_error=float(reference["overlap_reference_error"]))
        rows.append(row)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(args.output)


if __name__ == "__main__":
    main()
