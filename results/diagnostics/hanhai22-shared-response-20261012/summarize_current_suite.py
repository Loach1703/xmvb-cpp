#!/usr/bin/env python3
"""Summarize current NEO runs and saved same-input results, without reruns."""
import argparse
import csv
import hashlib
from pathlib import Path
import sys

from summarize_ciah_preconditioners import capture, read, row


def input_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else ""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--markdown", action="store_true")
    args = parser.parse_args()
    manifest = Path(__file__).with_name("preconditioner-suite.tsv")
    cases = list(csv.DictReader(manifest.open(), delimiter="\t"))
    results = []
    for case in cases:
        path = args.run_root / case["case"]
        result = row(case["case"], "hessian-diagonal", path)
        result["converged"] = result["converged"] and result["exit_code"] == "0"
        result["orbital_type"] = case["orbital_type"]
        result["input_sha256"] = input_hash(path / "input.xmi")
        result["physical_gradient"] = capture(
            r"Final physical \|g\|_\(G\^-1\)\s*:\s*(\S+)", read(path / "output.out"))
        result["baseline_steps"] = result["baseline_hvp"] = ""
        result["baseline_wall_s"] = result["baseline_energy_delta_Eh"] = ""
        result["same_input"] = ""
        if args.baseline:
            old_path = args.baseline / case["case"] / "hessian-diagonal"
            old = row(case["case"], "hessian-diagonal", old_path)
            old_hash = input_hash(old_path / "input.xmi") or input_hash(old_path.parent / "input.xmi")
            result["same_input"] = bool(old_hash and old_hash == result["input_sha256"])
            if result["same_input"] and old["converged"] and old["exit_code"] == "0":
                result["baseline_steps"] = old["steps"]
                result["baseline_hvp"] = old["recorded_orbital_hvp"]
                result["baseline_wall_s"] = old["wall_s"]
                if result["converged"] and result["energy_Eh"] and old["energy_Eh"]:
                    result["baseline_energy_delta_Eh"] = float(result["energy_Eh"]) - float(old["energy_Eh"])
        results.append(result)
    if not args.markdown:
        writer = csv.DictWriter(sys.stdout, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
        return
    print("| Case | Orbitals | Status | Steps | HVP | Energy / Eh | Wall / s | RSS / MiB | Old steps/HVP |")
    print("|---|---|---|---:|---:|---:|---:|---:|---|")
    for r in results:
        final = r["status"] == "finished"
        status = "converged" if r["converged"] else "failed/not converged" if final else r["status"]
        count = r["steps"] or (str(r["recorded_steps"]) + "+" if final else "—")
        hvp = r["recorded_orbital_hvp"] if r["converged"] else (
            str(r["recorded_orbital_hvp"]) + "+" if final else "—")
        wall = f'{r["wall_s"]:.2f}' if final and r["wall_s"] != "" else "—"
        rss = f'{r["peak_RSS_MiB"]:.1f}' if final and r["peak_RSS_MiB"] != "" else "—"
        old = f'{r["baseline_steps"]}/{r["baseline_hvp"]}' if r["baseline_steps"] else "—"
        print(f'| {r["system"]} | {r["orbital_type"]} | {status} | {count} | {hvp} | '
              f'{r["energy_Eh"] or "—"} | {wall} | {rss} | {old} |')
    print("\nOld results require identical input hashes and confirmed convergence; they are not rerun.")
    print("Nodes/builds and exact-exchange workspace budgets may differ: old wall times are contextual, not matched speedups.")
    print("A '+' marks accepted-step/action subtotals on failed runs, not complete attempt counts.")


if __name__ == "__main__":
    main()
