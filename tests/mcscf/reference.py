"""Export small independent CASCI references; run only inside Slurm."""
import os
from pathlib import Path

if not os.environ.get("SLURM_JOB_ID"):
    raise SystemExit("Run this molecular calculation in a Slurm compute job")

import numpy as np
from pyscf import ao2mo, gto, lib, mcscf, scf

lib.num_threads(int(os.environ.get("OMP_NUM_THREADS", "1")))

for name, atom, basis, electrons, orbitals, spin in [
    ("F2", "F 0 0 0; F 0 0 1.4", "sto-3g", 2, 2, 0),
    ("O2", "O 0 0 0; O 0 0 1.2", "sto-3g", 8, 6, 2),
    ("H4", "H 0 0 0; H 0 0 1.2; H 0 0 2.4; H 0 0 3.6", "sto-3g", 4, 4, 0),
]:
    mol = gto.M(atom=atom, basis=basis, spin=spin, verbose=0)
    hf = scf.ROHF(mol).run(conv_tol=1e-11)
    if not hf.converged:
        raise RuntimeError("reference HF did not converge")
    model = mcscf.CASCI(hf, orbitals, electrons)
    model.fcisolver.spin = spin
    model.fix_spin_(ss=spin/2*(spin/2+1))
    model.fcisolver.conv_tol = 1e-12
    energy = model.kernel()[0]
    h, core = model.get_h1eff(hf.mo_coeff)
    eri = ao2mo.restore(4, model.get_h2eff(hf.mo_coeff), orbitals)
    path = Path(name + ".integrals")
    with path.open("w") as stream:
        stream.write(f"{orbitals} {electrons} {spin} {core:.17g} {energy:.17g}\n")
        np.savetxt(stream, h, fmt="%.17g")
        np.savetxt(stream, eri, fmt="%.17g")
    print(name, "energy", energy, "spin", model.fcisolver.spin_square(
        model.ci, orbitals, model.nelecas)[0], flush=True)
