"""Check that fixed-orbital quintet CASCI contains the reference determinant."""

import json
import os
from pathlib import Path

assert os.environ.get("SLURM_JOB_ID"), "Run on a Slurm compute node"
import numpy as np
from pyscf import lib, mcscf, scf

lib.num_threads(int(os.environ["SLURM_CPUS_PER_TASK"]))
mol, saved = scf.chkfile.load_scf("quintet_rohf.chk")
mf = scf.ROHF(mol)
mf.mo_coeff, mf.mo_occ = saved["mo_coeff"], saved["mo_occ"]
mf.mo_energy, mf.e_tot = saved["mo_energy"], saved["e_tot"]
mf.converged = True
solver = mcscf.CASCI(mf, 11, (6, 2))
solver.max_memory = 120000
solver.fcisolver.conv_tol = 1e-10
solver.fcisolver.max_cycle = 200
solver.kernel(np.load("initial_mo.npy"))
assert solver.fcisolver.converged, "Fixed-orbital CI failed to converge"
assert solver.e_tot <= mf.e_tot + 1e-7, "Initial CAS lost the ROHF determinant"
report = dict(rohf_energy_Eh=float(mf.e_tot),
              initial_casci_energy_Eh=float(solver.e_tot),
              casci_minus_rohf_Eh=float(solver.e_tot - mf.e_tot))
Path("casci-check.json").write_text(json.dumps(report, indent=2) + "\n")
print("CASCI_CHECK", json.dumps(report), flush=True)
