#!/usr/bin/env python3
"""PySCF CASSCF(10,8) comparison for the native FeCl2 OEO input."""

from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import numpy as np
from pyscf import gto, mcscf, scf


def read_xmvb_orbitals(
    path: Path, n_ao: int, n_orbital: int, ao_scale: np.ndarray
) -> np.ndarray:
    text = path.read_text()
    block = text[text.lower().index("$gus"):].splitlines()[1:]
    header = re.compile(r"#\s*ORBITAL\s+(\d+)\s+NAO\s*=\s*(\d+)", re.I)
    pair = re.compile(
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?)\s+(\d+)"
    )
    orbitals: dict[int, np.ndarray] = {}
    current: int | None = None
    for line in block:
        match = header.search(line)
        if match:
            current = int(match.group(1)) - 1
            if int(match.group(2)) != n_ao:
                raise ValueError("unexpected AO dimension in $GUS")
            orbitals[current] = np.zeros(n_ao)
            continue
        if current is None:
            continue
        for value, index in pair.findall(line):
            orbitals[current][int(index) - 1] = float(value.replace("D", "E"))
    if sorted(orbitals) != list(range(n_orbital)):
        raise ValueError(f"expected {n_orbital} orbitals, found {sorted(orbitals)}")
    # xmvb-cpp uses individually normalized Cartesian AOs; PySCF/libcint does not.
    return np.column_stack([orbitals[i] for i in range(n_orbital)]) / ao_scale[:, None]


def complete_orthonormal_mos(
    overlap: np.ndarray, guess: np.ndarray, n_core: int
) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(overlap)
    sqrt_overlap = (eigenvectors * np.sqrt(eigenvalues)) @ eigenvectors.T
    inv_sqrt_overlap = (eigenvectors / np.sqrt(eigenvalues)) @ eigenvectors.T
    whitened = sqrt_overlap @ guess
    core, core_factor = np.linalg.qr(whitened[:, :n_core], mode="reduced")
    active_trial = whitened[:, n_core:] - core @ (core.T @ whitened[:, n_core:])
    active, active_factor = np.linalg.qr(active_trial, mode="reduced")
    print(f"CORE_MIN_QR_DIAG = {np.min(np.abs(np.diag(core_factor))):.8e}")
    print(f"ACTIVE_MIN_QR_DIAG = {np.min(np.abs(np.diag(active_factor))):.8e}")
    occupied = np.column_stack((core, active))
    complete, _ = np.linalg.qr(occupied, mode="complete")
    mo = inv_sqrt_overlap @ complete
    print(
        "MO_ORTHONORMALITY_ERROR = "
        f"{np.linalg.norm(mo.T @ overlap @ mo - np.eye(overlap.shape[0])):.8e}"
    )
    return mo


def write_xmvb_orbitals(
    source: Path, destination: Path, orbitals: np.ndarray
) -> None:
    """Replace the full-AO $GUS block using xmvb-cpp AO normalization."""
    text = source.read_text()
    start = text.lower().index("$gus")
    lines = [text[:start].rstrip(), "", "$GUS"]
    n_ao, n_orbital = orbitals.shape
    lines.append(" ".join([str(n_ao)] * n_orbital))
    for orbital in range(n_orbital):
        lines.append(f"# ORBITAL {orbital + 1:10d}  NAO = {n_ao:6d}")
        entries = [
            f"{orbitals[basis, orbital]: .16e} {basis + 1:6d}"
            for basis in range(n_ao)
        ]
        for first in range(0, n_ao, 4):
            lines.append(" ".join(entries[first:first + 4]))
    lines.append("$END")
    destination.write_text("\n".join(lines) + "\n")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: fecl2_oeo_pyscf.py FeCl2.xmi")
    mol = gto.M(
        atom=[
            ("Fe", (0.0, 0.0, 0.0)),
            ("Cl", (0.0, 0.0, 2.13631900)),
            ("Cl", (0.0, 0.0, -2.13631900)),
        ],
        basis="def2-tzvp",
        unit="Angstrom",
        charge=0,
        spin=4,
        cart=True,
        symmetry=False,
        verbose=4,
    )
    print(f"REFERENCE_NAO = {mol.nao_nr()}")
    overlap = mol.intor_symmetric("int1e_ovlp")
    ao_scale = np.sqrt(np.diag(overlap))
    guess = read_xmvb_orbitals(Path(sys.argv[1]), mol.nao_nr(), 33, ao_scale)
    mo = complete_orthonormal_mos(overlap, guess, 25)

    mean_field = scf.ROHF(mol)
    mean_field.max_memory = 24000
    casci = mcscf.CASCI(mean_field, 8, (7, 3))
    casci.fcisolver.spin = 4
    initial_energy = casci.kernel(mo)[0]
    print(f"REFERENCE_INITIAL_CASCI = {initial_energy:.15f}")

    history: list[tuple[int, float, float]] = []

    def callback(envs: dict[str, object]) -> None:
        macro = int(envs.get("imacro", len(history) + 1))
        energy = float(envs.get("e_tot", np.nan))
        gradient = float(envs.get("norm_gorb", np.nan))
        history.append((macro, energy, gradient))
        print(
            f"REFERENCE_MACRO = {macro} ENERGY = {energy:.15f} "
            f"GRAD = {gradient:.8e}"
        )

    casscf = mcscf.CASSCF(mean_field, 8, (7, 3))
    casscf.fcisolver.spin = 4
    casscf.conv_tol = 1.0e-7
    casscf.conv_tol_grad = 1.0e-3
    casscf.max_cycle_macro = 100
    casscf.callback = callback
    start = time.perf_counter()
    final_energy = casscf.kernel(mo)[0]
    final_xmvb_orbitals = casscf.mo_coeff[:, :33] * ao_scale[:, None]
    np.save("pyscf_final_mo.npy", casscf.mo_coeff)
    write_xmvb_orbitals(
        Path(sys.argv[1]), Path("FeCl2_pyscf.xmi"), final_xmvb_orbitals
    )
    print(f"REFERENCE_CASSCF = {final_energy:.15f}")
    print(f"REFERENCE_CONVERGED = {casscf.converged}")
    print(f"REFERENCE_MACRO_COUNT = {history[-1][0] if history else 0}")
    print(f"REFERENCE_WALL_SECONDS = {time.perf_counter() - start:.6f}")


if __name__ == "__main__":
    main()
