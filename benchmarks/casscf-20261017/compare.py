"""Summarize shared-guess energy parity without conflating solver iteration counts."""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("pyscf_result", type=Path)
    parser.add_argument("cpp_stdout", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    meta = json.loads(args.metadata.read_text())
    py = json.loads(args.pyscf_result.read_text())
    cpp = {}
    for line in args.cpp_stdout.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            cpp[key] = value
    required = ("initial_energy", "energy", "converged", "macros", "hvp_actions")
    if any(key not in cpp for key in required):
        raise ValueError("C++ runner did not emit all required fields")
    initial_gap = abs(float(cpp["initial_energy"]) - meta["initial_casci_energy"])
    final_gap = abs(float(cpp["energy"]) - py["energy"])
    result = dict(case=meta["case"], fixture_sha256=meta["fixture_sha256"],
                  py_solver=py.get("solver"), py_solver_mode=py.get("solver_mode"),
                  initial_casci_pyscf=meta["initial_casci_energy"],
                  initial_energy_cpp=float(cpp["initial_energy"]),
                  initial_absolute_error=initial_gap,
                  py_converged=py["converged"], cpp_converged=cpp["converged"] == "1",
                  py_energy=py["energy"], cpp_energy=float(cpp["energy"]),
                  final_absolute_error=final_gap,
                  py_macro_callbacks=py["macro_callback_count"],
                  cpp_macros=int(cpp["macros"]),
                  py_hvp_actions=py["hvp_actions"],
                  py_hvp_count_method=py.get("hvp_count_method"),
                  cpp_hvp_actions=int(cpp["hvp_actions"]),
                  initial_parity_1e_8=initial_gap <= 1e-8,
                  final_parity_1e_6=final_gap <= 1e-6)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if not result["initial_parity_1e_8"]:
        raise SystemExit("initial same-guess CASCI energy differs")
    if not result["py_converged"] or not result["cpp_converged"]:
        raise SystemExit("one or both CASSCF calculations did not converge")
    if not result["final_parity_1e_6"]:
        raise SystemExit("converged CASSCF energies differ")


if __name__ == "__main__":
    main()
