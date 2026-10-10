#!/usr/bin/env python3
"""Check distributed ERI fingerprints and five additive operator probes."""
import math
import pathlib
import re
import sys


def read(path):
    text = path.read_text()
    return dict(re.findall(r"\b([a-z_0-9]+)=([^\s]+)", text))


root = pathlib.Path(sys.argv[1])
reference = read(root / "1x32/output.out")
for layout in ("2x16", "2x32"):
    measured = read(root / layout / "output.out")
    for key in ("total_integrals", "hash_sum", "hash_xor"):
        if measured[key] != reference[key]:
            raise RuntimeError(f"{layout}: ERI ownership/fingerprint mismatch: {key}")
    for channel in range(5):
        key = f"probe_{channel}"
        if not math.isclose(float(measured[key]), float(reference[key]),
                            rel_tol=2e-12, abs_tol=1e-8):
            raise RuntimeError(f"{layout}: distributed contraction mismatch: {key}")
    print(f"{layout}: exact integral count/fingerprints and all operator probes passed")
