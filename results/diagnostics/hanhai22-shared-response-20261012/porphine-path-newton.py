"""One-macro coupled PySCF control for the original singlet guess."""
import json
import os
from pathlib import Path
import runpy
import sys

assert os.environ.get("SLURM_JOB_ID"), "Use a Slurm compute node"
import numpy as np
from pyscf.mcscf import newton_casscf

root = Path(sys.argv[1])
original_kernel = newton_casscf.CASSCF.kernel


def traced_kernel(self, *args, **kwargs):
    self.max_cycle_macro = 1
    np.save("initial-mo.npy", args[0])
    previous_callback = self.callback

    def snapshot(env):
        # newton_casscf's callback casdm1 is stale until post-loop finalization.
        # Contract the accepted CI vector; do not use that stale density.
        density = self.fcisolver.make_rdm1(env["fcivec"], self.ncas, self.nelecas)
        np.savez("macro-001.npz", mo=env["mo"], active_density=density,
                 macro=env["imacro"], energy=env["e_tot"])
        record = dict(energy=float(env["e_tot"]), macro=int(env["imacro"]),
                      coupled_gradient_norm=float(env["norm_gall"]),
                      hvp=int(env["tot_hop"]), keyframes=int(env["tot_kf"]))
        Path("path.json").write_text(json.dumps(record, indent=2)+"\n")
        if callable(previous_callback):
            previous_callback(env)

    self.callback = snapshot
    try:
        return original_kernel(self, *args, **kwargs)
    finally:
        self.callback = previous_callback


newton_casscf.CASSCF.kernel = traced_kernel
sys.argv = [str(root / "porphine-cross-restart.py"), str(root),
            "newton", "converge", "original", "singlet"]
runpy.run_path(str(root / "porphine-cross-restart.py"), run_name="__main__")
