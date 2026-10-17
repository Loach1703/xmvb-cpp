"""Recover the established FeCl2 full-MO seed from its $GUS and Molden AO layout.

Adapted from results/diagnostics/fecl2-pyscf-reference-20261011/newton/
FeCl2_pyscf.py. The first 33 orbitals come from the XMVB deck; the remaining
virtual complement is completed deterministically in the AO-overlap metric.
"""

from pathlib import Path
import re

import numpy as np


def read_xmvb_orbitals(path: Path, nbf: int, count: int, ao_scale):
    contents = path.read_text()
    block = contents[contents.lower().index("$gus"):].splitlines()[1:]
    header = re.compile(r"#\s*ORBITAL\s+(\d+)\s+NAO\s*=\s*(\d+)", re.I)
    pair = re.compile(r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?)\s+(\d+)")
    orbitals = {}
    current = None
    for line in block:
        match = header.search(line)
        if match:
            current = int(match.group(1)) - 1
            if int(match.group(2)) != nbf:
                raise ValueError("$GUS AO dimension differs from Molden layout")
            orbitals[current] = np.zeros(nbf)
            continue
        if current is not None:
            for value, index in pair.findall(line):
                orbitals[current][int(index) - 1] = float(value.replace("D", "E"))
    if sorted(orbitals) != list(range(count)):
        raise ValueError("$GUS does not contain the expected core and active columns")
    # XMVB stores individually normalized Cartesian AO coefficients, while
    # the PySCF/Molden AO basis uses the unscaled Cartesian functions.
    return np.column_stack([orbitals[i] for i in range(count)]) / ao_scale[:, None]


def complete_orthonormal_mos(overlap, guess, ncore):
    eigenvalues, eigenvectors = np.linalg.eigh(overlap)
    if eigenvalues[0] <= 0:
        raise ValueError("AO overlap is not positive definite")
    sqrt_overlap = (eigenvectors * np.sqrt(eigenvalues)) @ eigenvectors.T
    inv_sqrt_overlap = (eigenvectors / np.sqrt(eigenvalues)) @ eigenvectors.T
    whitened = sqrt_overlap @ guess
    core, core_factor = np.linalg.qr(whitened[:, :ncore], mode="reduced")
    active_trial = whitened[:, ncore:] - core @ (core.T @ whitened[:, ncore:])
    active, active_factor = np.linalg.qr(active_trial, mode="reduced")
    if min(np.min(np.abs(np.diag(core_factor))),
           np.min(np.abs(np.diag(active_factor)))) < 1e-8:
        raise ValueError("archived FeCl2 core/active guess is rank deficient")
    complete, _ = np.linalg.qr(np.column_stack((core, active)), mode="complete")
    mo = inv_sqrt_overlap @ complete
    if np.max(np.abs(mo.T @ overlap @ mo - np.eye(overlap.shape[0]))) > 1e-8:
        raise ValueError("completed FeCl2 MO frame is not S-orthonormal")
    return mo
