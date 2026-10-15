"""Default-AH CASSCF from export_shared.py artifacts; no new HF initial guess."""
import hashlib
import json
import os
from pathlib import Path
import time


def main():
    assert os.environ.get("SLURM_JOB_ID"), "Run on a compute node"
    import numpy as np
    from pyscf import lib, scf, mcscf, fci
    lib.num_threads(int(os.environ["SLURM_CPUS_PER_TASK"]))
    data = json.loads(Path("shared.json").read_text())
    for name, digest in data["sha256"].items():
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == digest
    mol, _ = scf.chkfile.load_scf("guess.chk")
    mo = np.load("initial_mo.npy",allow_pickle=False)
    assert mol.cart and not mol.symmetry and not mol.has_ecp() and mol.spin == data["spin"]
    spin = data["spin"]
    ne = data["nelecas"]
    mc = mcscf.CASSCF(scf.ROHF(mol),data["ncas"],((ne+spin)//2,(ne-spin)//2))
    mc.fix_spin_(ss=(spin/2)*(spin/2+1))
    mc.conv_tol, mc.conv_tol_grad = 1e-7, 1e-3
    mc.fcisolver.conv_tol, mc.fcisolver.max_cycle = 1e-10, 200
    mc.max_cycle_macro, mc.max_memory = 1000, 120000
    mc.chkfile, mc.verbose = "optimized.chk", 5
    start = time.monotonic()
    mc.kernel(mo)
    result = dict(converged=bool(mc.converged), energy=float(mc.e_tot),
                  elapsed_seconds=time.monotonic()-start,
                  spin_square=float(fci.spin_op.spin_square(mc.ci,mc.ncas,mc.nelecas)[0]))
    Path("result.json").write_text(json.dumps(result,indent=2)+"\n")
    np.save("final_mo.npy",mc.mo_coeff)


if __name__ == "__main__":
    main()
