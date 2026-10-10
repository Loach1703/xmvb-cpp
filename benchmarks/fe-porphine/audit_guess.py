"""Audit the saved initial frame without rerunning SCF or orbital optimization."""

import json
import os
from pathlib import Path
import sys

assert os.environ.get("SLURM_JOB_ID"), "Run on a Slurm compute node"
import numpy as np
from pyscf import lib, scf

lib.num_threads(int(os.environ["SLURM_CPUS_PER_TASK"]))
prepared = Path(sys.argv[1])
mol, saved = scf.chkfile.load_scf(str(prepared / "quintet_rohf.chk"))
mo = np.load(prepared / "initial_mo.npy")
reference = saved["mo_coeff"]
occupations = saved["mo_occ"]
s = mol.intor_symmetric("int1e_ovlp")
n_core = (mol.nelectron - 8) // 2
core, active = mo[:, :n_core], mo[:, n_core:n_core + 11]
occupied = reference[:, occupations > 0]
closed = reference[:, occupations == 2]
singly = reference[:, occupations == 1]
carrier = mo[:, :n_core + 11]
carrier_weights = carrier.T @ s @ occupied
core_in_closed = closed.T @ s @ core
core_in_singly = singly.T @ s @ core
report = {
    "n_bf": mol.nao_nr(),
    "n_core": n_core,
    "rohf_energy_Eh": float(saved["e_tot"]),
    "initial_core_population": float(np.sum(
        (core.T @ s @ reference) ** 2 * occupations[None, :])),
    "initial_active_population": float(np.sum(
        (active.T @ s @ reference) ** 2 * occupations[None, :])),
    "occupied_leakage_outside_core_and_active": float(
        occupied.shape[1] - np.sum(carrier_weights ** 2)),
    "core_leakage_outside_rohf_closed_span": float(
        n_core - np.sum(core_in_closed ** 2)),
    "rohf_singly_occupied_weight_in_core": float(np.sum(core_in_singly ** 2)),
    "max_orbital_gram_error": float(np.max(np.abs(mo.T @ s @ mo - np.eye(mo.shape[1])))),
}
print(json.dumps(report, indent=2))
