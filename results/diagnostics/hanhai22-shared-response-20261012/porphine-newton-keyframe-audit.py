"""One-macro PySCF coupled-Newton audit of the first keyframe gradient.

Run through the same Slurm environment and original singlet input as
porphine-path-newton.py.  The extra exact gradient uses the keyframe CI vector
without a CASCI solve and is never returned to the Newton optimizer.
"""

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
original_extract_rotation = newton_casscf.extract_rotation
original_gen_g_hop = newton_casscf.gen_g_hop
audit = {"pending": False, "recorded": False, "gen_calls": 0}


def audited_extract_rotation(casscf, dr, u, ci0):
    rotated_u, rotated_ci = original_extract_rotation(casscf, dr, u, ci0)
    if not audit["recorded"] and not audit["pending"]:
        audit["dr"] = np.array(dr, copy=True)
        audit["u"] = np.array(rotated_u, copy=True)
        audit["ci"] = np.array(rotated_ci, copy=True)
        audit["pending"] = True
    return rotated_u, rotated_ci


def audited_gen_g_hop(casscf, mo, ci0, eris, verbose=None):
    audit["gen_calls"] += 1
    assert audit["gen_calls"] == 1, "Expected exactly one Newton macro"
    g0, g_update, h_op, h_diag = original_gen_g_hop(
        casscf, mo, ci0, eris, verbose=verbose)
    initial_gradient = np.array(g0, copy=True)
    initial_mo = np.array(mo, copy=True)
    n_orbital = int(np.count_nonzero(casscf.uniq_var_indices(
        mo.shape[1], casscf.ncore, casscf.ncas, casscf.frozen)))

    def audited_g_update(u, ci):
        approximate = g_update(u, ci)
        if audit["pending"]:
            assert not audit["recorded"]
            assert np.array_equal(u, audit["u"])
            assert np.array_equal(ci, audit["ci"])
            dr = audit["dr"]
            predicted = initial_gradient + h_op(dr)
            keyframe_mo = initial_mo @ u
            keyframe_eris = casscf.ao2mo(keyframe_mo)
            # Call the saved function directly: the diagnostic must not recurse
            # through this module's gen_g_hop hook or reoptimize the CI vector.
            exact = original_gen_g_hop(
                casscf, keyframe_mo, ci, keyframe_eris)[0]
            assert predicted.shape == approximate.shape == exact.shape

            def split_norms(vector):
                return dict(total=float(np.linalg.norm(vector)),
                            orbital=float(np.linalg.norm(vector[:n_orbital])),
                            ci=float(np.linalg.norm(vector[n_orbital:])))

            record = dict(
                macro=1, keyframe=1, n_orbital=n_orbital,
                n_ci=int(exact.size - n_orbital),
                dr=split_norms(dr),
                predicted=split_norms(predicted),
                approximate=split_norms(approximate),
                exact=split_norms(exact),
                predicted_minus_approximate=split_norms(predicted - approximate),
                predicted_minus_exact=split_norms(predicted - exact),
                approximate_minus_exact=split_norms(approximate - exact),
            )
            np.savez("first-keyframe.npz", dr=dr, u=u, ci=ci)
            Path("first-keyframe.json").write_text(json.dumps(record, indent=2) + "\n")
            print("FIRST_KEYFRAME_AUDIT", json.dumps(record), flush=True)
            audit["pending"] = False
            audit["recorded"] = True
        return approximate

    return g0, audited_g_update, h_op, h_diag


def traced_kernel(self, *args, **kwargs):
    self.max_cycle_macro = 1
    np.save("initial-mo.npy", args[0])
    previous_callback = self.callback

    def snapshot(env):
        # newton_casscf's callback casdm1 is stale until post-loop finalization.
        density = self.fcisolver.make_rdm1(env["fcivec"], self.ncas, self.nelecas)
        np.savez("macro-001.npz", mo=env["mo"], active_density=density,
                 macro=env["imacro"], energy=env["e_tot"])
        record = dict(energy=float(env["e_tot"]), macro=int(env["imacro"]),
                      coupled_gradient_norm=float(env["norm_gall"]),
                      hvp=int(env["tot_hop"]), keyframes=int(env["tot_kf"]))
        Path("path.json").write_text(json.dumps(record, indent=2) + "\n")
        if callable(previous_callback):
            previous_callback(env)

    self.callback = snapshot
    newton_casscf.extract_rotation = audited_extract_rotation
    newton_casscf.gen_g_hop = audited_gen_g_hop
    try:
        result = original_kernel(self, *args, **kwargs)
        assert audit["recorded"], "No intermediate keyframe was audited"
        return result
    finally:
        newton_casscf.gen_g_hop = original_gen_g_hop
        newton_casscf.extract_rotation = original_extract_rotation
        self.callback = previous_callback


newton_casscf.CASSCF.kernel = traced_kernel
sys.argv = [str(root / "porphine-cross-restart.py"), str(root),
            "newton", "converge", "original", "singlet"]
runpy.run_path(str(root / "porphine-cross-restart.py"), run_name="__main__")
