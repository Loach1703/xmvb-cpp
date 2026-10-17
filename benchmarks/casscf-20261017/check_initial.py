"""Gate a full CASSCF run on exact same-frame initial CASCI parity."""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("cpp_stdout", type=Path)
    args = parser.parse_args()
    meta = json.loads(args.metadata.read_text())
    lines = args.cpp_stdout.read_text().splitlines()
    energies = [float(line.split("=", 1)[1]) for line in lines
                if line.startswith("initial_energy=")]
    if len(energies) != 1:
        raise ValueError("expected exactly one C++ initial CASCI energy")
    difference = abs(energies[0] - meta["initial_casci_energy"])
    print(f"initial CASCI |C++ - PySCF| = {difference:.12g} Eh", flush=True)
    if difference > 1e-8:
        raise SystemExit("initial CASCI parity failed; skip optimization")


if __name__ == "__main__":
    main()
