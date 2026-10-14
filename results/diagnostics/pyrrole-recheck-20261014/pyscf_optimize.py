"""Two-state SA-CASSCF optimization from the new attached XMVB $GUS."""
import hashlib
import json
import os
import re
import time
from pathlib import Path

import numpy as np
import pyscf
from pyscf import fci, lib, mcscf, scf
from pyscf.tools import molden

assert os.environ.get("SLURM_JOB_ID"), "Run on a Slurm compute node"
lib.num_threads(int(os.environ["SLURM_CPUS_PER_TASK"]))
source = Path("input.xmi").read_text()
print("INPUT_SHA256", hashlib.sha256(source.encode()).hexdigest(), flush=True)
print("PYSCF_VERSION", pyscf.__version__, flush=True)

# Use only the basis/shell layout from Molden, NOT its final orbital columns.
# This layout was validated in the preceding fixed-orbital cross-check.
mol = molden.load("basis_layout.molden")[0]
mol.charge = 1
mol.spin = 1
mol.symmetry = False
mol.verbose = 4
mol.max_memory = 48000
assert mol.cart and mol.nao_nr() == 340 and mol.nelectron == 71
overlap = mol.intor_symmetric("int1e_ovlp")
ao_scale = np.sqrt(overlap.diagonal())

# $GUS coefficients use individually normalized Cartesian AOs. PySCF uses
# radial-normalized Cartesian functions with the same shell-local xx,xy,xz,
# yy,yz,zz ordering. Read all coefficients at their original input precision.
guess = np.zeros((340, 40))
seen = set()
orbital = None
header = re.compile(r"#\s*ORBITAL\s+(\d+)\s+NAO\s*=\s*(\d+)", re.I)
pair = re.compile(r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?)\s+(\d+)")
for line in re.split(r"\$GUS", source, flags=re.I)[1].splitlines():
    if line.strip().upper() == "$END":
        break
    match = header.search(line)
    if match:
        orbital = int(match[1]) - 1
        assert int(match[2]) == 340 and 0 <= orbital < 40
        seen.add(orbital)
    elif orbital is not None:
        for value, index in pair.findall(line):
            guess[int(index) - 1, orbital] = float(value.replace("D", "E"))
assert seen == set(range(40))
guess /= ao_scale[:, None]

# Orthogonal representatives preserve the occupied core span and the active
# span modulo that core. This is a gauge conversion, not a new SCF initial guess.
values, vectors = np.linalg.eigh(overlap)
sqrt_s = (vectors * np.sqrt(values)) @ vectors.T
inv_sqrt_s = (vectors / np.sqrt(values)) @ vectors.T
white = sqrt_s @ guess
core, core_factor = np.linalg.qr(white[:, :30], mode="reduced")
active_trial = white[:, 30:40] - core @ (core.T @ white[:, 30:40])
active, active_factor = np.linalg.qr(active_trial, mode="reduced")
print("CORE_MIN_QR_DIAG", np.min(np.abs(np.diag(core_factor))), flush=True)
print("ACTIVE_MIN_QR_DIAG", np.min(np.abs(np.diag(active_factor))), flush=True)
assert np.min(np.abs(np.diag(active_factor))) > 1e-8
full = np.linalg.qr(np.column_stack((core, active)), mode="complete")[0]
mo = inv_sqrt_s @ full
error = np.max(np.abs(mo.T @ overlap @ mo - np.eye(340)))
print("INITIAL_ORTHONORMALITY_MAX", error, flush=True)
assert error < 1e-8
np.save("initial_mo.npy", mo)

def doublet_solver():
    solver = fci.addons.fix_spin_(fci.direct_spin1.FCISolver(mol), shift=1.0, ss=0.75)
    solver.conv_tol = 1e-11
    solver.max_cycle = 300
    return solver

mean_field = scf.ROHF(mol)  # No SCF kernel: never replace the supplied initial guess.
initial = mcscf.CASCI(mean_field, 10, (6, 5))
initial.canonicalization = False
initial.fcisolver = doublet_solver()
initial.fcisolver.nroots = 2
initial_energies = np.asarray(initial.kernel(mo)[0])
average = float(initial_energies.mean())
print("INITIAL_ROOT_ENERGIES", initial_energies.tolist(), flush=True)
print("INITIAL_AVERAGE_ENERGY", average, flush=True)
print("CPP_INITIAL_AVERAGE_ENERGY", -417.496351311941, flush=True)
assert abs(average - (-417.496351311941)) < 2e-6, "Initial energy mismatch: do not optimize"

optimizer = mcscf.CASSCF(mean_field, 10, (6, 5))
optimizer.fcisolver = doublet_solver()
optimizer = optimizer.state_average_([0.5, 0.5])
optimizer.max_cycle_macro = 200
optimizer.conv_tol = 1e-8
optimizer.conv_tol_grad = 1e-5
optimizer.chkfile = "pyscf.chk"
history = []

def callback(env):
    row = {"macro": int(env.get("imacro", -1)),
           "energy": float(env.get("e_tot", np.nan)),
           "gradient": float(env.get("norm_gorb", np.nan)),
           "elapsed_s": time.perf_counter() - start}
    history.append(row)
    print("MACRO_TRACE", json.dumps(row), flush=True)
    Path("history.json").write_text(json.dumps(history, indent=2))

optimizer.callback = callback
start = time.perf_counter()
optimizer.kernel(mo, ci0=initial.ci)
result = {"converged": bool(optimizer.converged),
          "initial_roots": initial_energies.tolist(),
          "initial_average": average,
          "final_roots": np.asarray(optimizer.e_states).tolist(),
          "final_average": float(optimizer.e_tot),
          "spin_squared": [float(fci.spin_op.spin_square(ci, 10, (6, 5))[0]) for ci in optimizer.ci],
          "optimization_wall_s": time.perf_counter() - start,
          "conv_tol": optimizer.conv_tol, "conv_tol_grad": optimizer.conv_tol_grad,
          "history": history}
np.save("final_mo.npy", optimizer.mo_coeff)
molden.from_mo(mol, "final.molden", optimizer.mo_coeff)
Path("result.json").write_text(json.dumps(result, indent=2))
print("FINAL_RESULT", json.dumps(result), flush=True)
