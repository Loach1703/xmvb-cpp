"""Recover raw XMVB coefficients from its Molden export, without AO algebra."""

import os
from pathlib import Path
import re
import sys

assert os.environ.get("SLURM_JOB_ID"), "Run the diagnostic under Slurm"
original, molden, target = map(Path, sys.argv[1:])
text = molden.read_text()
gto = re.search(r"\[GTO\](.*?)(?=\[)", text, re.S | re.I)[1]
# Inverse of output/molden/writer.cpp's shell-local permutation. The export
# writes stored raw coefficients directly; do not introduce an AO rescaling.
orders = {"s": [0], "p": [0, 1, 2], "d": [0, 3, 5, 1, 2, 4],
          "f": [0, 6, 9, 3, 1, 2, 5, 8, 7, 4]}
permutation = []
for shell in re.finditer(r"^[ \t]*([spdf])[ \t]+\d+(?:[ \t]+[\d.Ee+-]+)?[ \t]*$", gto, re.M | re.I):
    order = orders[shell[1].lower()]
    offset = len(permutation)
    permutation.extend(offset + i for i in order)
assert len(permutation) == 469, "Unexpected basis or shell convention"
mo = re.split(r"\[MO\]", text, flags=re.I)[1]
columns = []
current = None
for line in mo.splitlines():
    if re.match(r"\s*Sym\s*=", line, re.I):
        if current is not None:
            assert all(value is not None for value in current)
            columns.append(current)
        current = [None] * len(permutation)
    match = re.fullmatch(r"\s*(\d+)\s+([\d.Ee+-]+)\s*", line)
    if match and current is not None:
        index = int(match[1]) - 1
        current[permutation[index]] = float(match[2])
if current is not None:
    assert all(value is not None for value in current)
    columns.append(current)
assert len(columns) == 100, "Expected 89 inactive and 11 active orbitals"
prefix = re.split(r"\$GUS", original.read_text(), flags=re.I)[0].rstrip()
lines = [prefix, "$GUS", " ".join(["469"] * 100)]
for orbital, column in enumerate(columns, 1):
    lines.append(f"# ORBITAL {orbital} NAO = 469")
    for first in range(0, len(column), 4):
        lines.append(" ".join(f"{column[i]:.17e} {i+1}"
                              for i in range(first, min(first+4, len(column)))))
lines.append("$END")
target.write_text("\n".join(lines) + "\n")
print(f"Converted {len(columns)} raw orbitals, {len(permutation)} AO entries to {target}")
