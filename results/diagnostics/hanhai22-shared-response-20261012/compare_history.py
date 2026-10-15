#!/usr/bin/env python3
"""Pair current native NEO artifacts with saved c2776c4 results; no reruns."""
import argparse
import csv
from pathlib import Path
import sys

from summarize_ciah_preconditioners import capture, read, row


ALIASES = {
    "241": "241_vb", "240": "240_vb", "7975": "7975_vb", "7963": "7963_vb",
    "YAMSAI": "csd_YAMSAI_vb", "NEGPIA": "csd_NEGPIA_vb",
    "QOCGUM": "csd_QOCGUM_vb", "67970": "67970_vb",
    "CERRAS": "csd_CERRAS_vb", "LOFLEA": "csd_LOFLEA_vb",
    "10698": "10698_VBSCF", "F2-OEO": "F2_OEO",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--history", type=Path,
                        default=Path(__file__).parents[2] / "history/paired_runs.csv")
    parser.add_argument("--previous", type=Path)
    args = parser.parse_args()
    baseline = {r["case"]: r for r in csv.DictReader(args.history.open())
                if r["method"] == "neo" and r["revision"] == "c2776c4"}
    cases = list(csv.DictReader(Path(__file__).with_name("preconditioner-suite.tsv").open(),
                               delimiter="\t"))
    rows = []
    for case in cases:
        name = case["case"]
        r = row(name, "hessian-diagonal", args.run_root / name)
        old = baseline[ALIASES.get(name, name)]
        text = read(args.run_root / name / "output.out")
        total = capture(r"End-to-end wall time\s*:\s*([\d.]+)", text)
        initial = capture(r"Initial total energy\s*:\s*(-?[\d.]+)", text)
        converged = r["converged"] and r["exit_code"] == "0"
        old_total = float(old["total_wall_s"])
        energy_delta = float(r["energy_Eh"]) - float(old["energy_final_hartree"]) \
            if converged else ""
        previous_text = read(args.previous / name / "output.out") if args.previous else ""
        previous_total = capture(r"End-to-end wall time\s*:\s*([\d.]+)", previous_text)
        rows.append({
            "case": name, "status": "converged" if converged else
                "not_converged" if r["status"] == "finished" else r["status"],
            "exit_code": r["exit_code"], "termination": r["termination"],
            "old_steps": old["accepted_steps"], "steps": r["steps"],
            "old_hvp": old["orbital_hvp_actions"], "hvp": r["recorded_orbital_hvp"],
            "old_total_s": old_total, "total_s": float(total) if converged and total else "",
            "time_ratio": float(total)/old_total if converged and total else "",
            "previous_total_s": previous_total,
            "peak_RSS_MiB": r["peak_RSS_MiB"] if converged else "",
            "energy_Eh": r["energy_Eh"] if converged else "",
            "energy_delta_Eh": energy_delta,
            "initial_energy_delta_Eh": float(initial)-float(old["energy_initial_hartree"])
                if initial else "",
            "node": r["node"], "note": "historical FeCl2 endpoint differs from PySCF reference"
                if name == "FeCl2" else "historical Hanhai25 versus current Hanhai22",
        })
    writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
