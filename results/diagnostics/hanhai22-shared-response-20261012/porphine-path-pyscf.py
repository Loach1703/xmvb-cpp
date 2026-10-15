"""Record the existing original-guess AH path without changing its solver."""
import json
import os
from pathlib import Path
import runpy
import sys

assert os.environ.get("SLURM_JOB_ID"), "Use a Slurm compute node"
import numpy as np
from pyscf.mcscf import mc1step

root = Path(sys.argv[1])
original_kernel = mc1step.CASSCF.kernel


def traced_kernel(self, *args, **kwargs):
    """Add passive snapshots; the four-macro cap is diagnostic only."""
    self.max_cycle_macro = 4
    initial_mo = args[0] if args else kwargs.get("mo_coeff", self.mo_coeff)
    np.save("initial-mo.npy", initial_mo)
    previous_callback = self.callback

    def snapshot(env):
        macro, micro = int(env["imacro"]), int(env["imicro"])
        accepted = env.get("rota") is None
        tag = f"macro-{macro:03d}" if accepted else f"micro-{macro:03d}-{micro:03d}"
        mo = env["mo"] if accepted else env["mo"] @ env["u"]
        # The micro density is approximate. Only accepted macro densities
        # follow the fully solved CI problem; never label them interchangeably.
        np.savez(f"{tag}.npz", mo=mo, active_density=env["casdm1"],
                 macro=macro, micro=micro, accepted=accepted,
                 energy=env["e_tot"] if accepted else np.nan)
        record = dict(tag=tag, accepted=accepted,
                      energy=float(env["e_tot"]) if accepted else None,
                      orbital_gradient_norm=float(env["norm_gorb"]),
                      density_change=float(env["norm_ddm"]),
                      rotation_norm=float(env["norm_t"]))
        with open("path.jsonl", "a") as file:
            file.write(json.dumps(record) + "\n")
        if callable(previous_callback):
            previous_callback(env)

    self.callback = snapshot
    try:
        return original_kernel(self, *args, **kwargs)
    finally:
        self.callback = previous_callback


mc1step.CASSCF.kernel = traced_kernel
sys.argv = [str(root / "porphine-cross-restart.py"), str(root),
            "default-ah", "converge", "original", "singlet"]
runpy.run_path(str(root / "porphine-cross-restart.py"), run_name="__main__")
