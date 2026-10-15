#!/usr/bin/env python3
"""Report matched CIAH preconditioner runs without treating partial runs as final."""

import argparse
import csv
from pathlib import Path
import sys

from summarize_ciah_preconditioners import row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--markdown", action="store_true")
    args = parser.parse_args()
    manifest = Path(__file__).with_name("preconditioner-suite.tsv")
    cases = list(csv.DictReader(manifest.open(), delimiter="\t"))
    rows = []
    for case in cases:
        for preconditioner in ("hessian-diagonal", "one-electron"):
            result = row(case["case"], preconditioner,
                         args.run_root / case["case"] / preconditioner)
            result["orbital_type"] = case["orbital_type"]
            rows.append(result)
    if not args.markdown:
        writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
        return
    print("| System | Preconditioner | Status | Steps | HVP | Energy (Eh) | Wall (s) | RSS (MiB) |")
    print("|---|---|---|---:|---:|---:|---:|---:|")
    for result in rows:
        status = ("converged" if result["converged"] else
                  "not converged" if result["status"] == "finished" else result["status"])
        # Counters on an unfinished run are partial; do not report them as totals.
        final = result["status"] == "finished"
        hvp = result["recorded_orbital_hvp"] if final else "—"
        if final and not result["converged"] and hvp != "":
            hvp = f"{hvp}+"
        wall = f'{result["wall_s"]:.2f}' if final and result["wall_s"] != "" else "—"
        rss = f'{result["peak_RSS_MiB"]:.1f}' if final and result["peak_RSS_MiB"] != "" else "—"
        print(f'| {result["system"]} | {result["preconditioner"]} | {status} | '
              f'{result["steps"] or (result["recorded_steps"] if final else "") or "—"} | '
              f'{hvp} | {result["energy_Eh"] or "—"} | {wall} | {rss} |')
    print("\nOn failed runs, steps are accepted steps; HVP marked '+' is only the recorded subtotal.")
    print("\nPaired energy differences (one-electron minus analytic diagonal):")
    for analytic, one_electron in zip(rows[::2], rows[1::2]):
        if analytic["converged"] and one_electron["converged"]:
            delta = float(one_electron["energy_Eh"]) - float(analytic["energy_Eh"])
            print(f'- {analytic["system"]}: {delta:+.12g} Eh')


if __name__ == "__main__":
    main()
