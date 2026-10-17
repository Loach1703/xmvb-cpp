#!/usr/bin/env python3
"""Analyze frozen selected-CI columns without rerunning a physical operator.

Distance is half the L1 difference of two structures' spatial occupations,
not a determinant excitation rank. The inverse model is occupation-block
diagonal, so complete distance bands are orthogonal in its inverse norm.
Individual row products need not be positive; sum a complete band first.
The current binary capture contract is little-endian uint64/int32/float64.
"""

import argparse
from array import array
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import sys


def read_array(path, kind, size):
    result = array(kind)
    result.frombytes(path.read_bytes())
    if sys.byteorder != "little":
        result.byteswap()
    if len(result) != size:
        raise ValueError(f"unexpected payload size: {path}")
    return result


def analyze(root):
    lines = (root / "models.out").read_text().splitlines()
    n = int(next(line.split("=", 1)[1] for line in lines
                 if line.startswith("model_dimension=")))
    table = [line for line in lines if line.startswith("diagonal_sample,")]
    samples = []
    for row in csv.DictReader(io.StringIO("\n".join(table))):
        row.pop("diagonal_sample")
        samples.append({name: int(value) if name in ("structure", "ionic_order")
                        else float(value) for name, value in row.items()})
    occupations = read_array(root / "probe.occupations.bin", "Q", 2 * n)
    electrons = {occupations[2 * i].bit_count() + occupations[2 * i + 1].bit_count()
                 for i in range(n)}
    if len(electrons) != 1:
        raise ValueError("structure occupations do not conserve electron number")
    labels = read_array(root / "probe.structures.bin", "i", len(samples))
    if [i + 1 for i in labels] != [row["structure"] for row in samples]:
        raise ValueError("probe labels disagree with the scalar report")
    bands = []
    for sample, label in zip(samples, labels):
        prefix = root / f"probe.{label + 1}"
        shifted = read_array(Path(str(prefix) + ".shifted.bin"), "d", n)
        inverse = read_array(Path(str(prefix) + ".inverse.bin"), "d", n)
        occupied, paired = occupations[2 * label:2 * label + 2]
        products = {}
        for row in range(n):
            l1 = ((occupied ^ occupations[2 * row]).bit_count() +
                  (paired ^ occupations[2 * row + 1]).bit_count())
            if l1 % 2:
                raise ValueError("structure occupations do not conserve electron number")
            products.setdefault(l1 // 2, []).append(shifted[row] * inverse[row])
        norms = {distance: math.fsum(values) for distance, values in products.items()}
        total = math.fsum(norms.values())
        if not total > 0:
            raise ValueError("probe inverse norm is not positive")
        cumulative = 0.0
        for distance, squared in sorted(norms.items()):
            cumulative += squared / total
            bands.append(dict(structure=label + 1, ionic_order=sample["ionic_order"],
                distance=distance, rows=len(products[distance]),
                squared_inverse_norm_fraction=squared / total,
                cumulative_squared_inverse_norm_fraction=cumulative))
        sample["relative_gap_model_error"] = (
            sample["model_gap_per_G"] / sample["exact_gap_per_G"] - 1)
    hashes = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
              for path in sorted(root.glob("probe.*.bin"))}
    return dict(dimension=n, electrons=electrons.pop(), sampled_structures=len(samples),
                diagonal_samples=samples, distance_bands=bands, probe_sha256=hashes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze(args.root), indent=2))
