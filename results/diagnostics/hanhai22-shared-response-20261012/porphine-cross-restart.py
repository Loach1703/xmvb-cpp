"""Two-macrostep PySCF restart from the actual CPP triplet wavefunction."""
import json
import os
from pathlib import Path
import sys
import time

assert os.environ.get("SLURM_JOB_ID"), "Run on a Slurm compute node"
import numpy as np
from scipy.linalg import null_space, solve_triangular
from pyscf import fci, lib, mcscf, scf

lib.num_threads(int(os.environ["SLURM_CPUS_PER_TASK"]))
root, mode = Path(sys.argv[1]), sys.argv[2]
source = root / "porphine-coupled-audit-1247963"
prefix = source / "audit_orbital_gradient_"
mol, saved = scf.chkfile.load_scf(str(root / "porphine-reference.chk"))
mol.spin = 2
mol.build(False, False)
assert mol.cart and mol.nao_nr() == 469 and mol.nelectron == 186 and not mol._ecp
mf = scf.ROHF(mol)
mf.mo_coeff, mf.mo_occ = saved["mo_coeff"], saved["mo_occ"]
mf.mo_energy, mf.e_tot, mf.converged = saved["mo_energy"], saved["e_tot"], True
s = mf.get_ovlp()
d = np.sqrt(s.diagonal())
cpp_s = np.fromfile(str(prefix) + "ao_overlap_f64.bin").reshape((469, 469), order="F")
ao_error = float(np.max(np.abs(cpp_s - s / d[:, None] / d[None, :])))
assert ao_error < 1e-8, f"AO normalization/order mismatch: {ao_error}"
raw = np.fromfile(str(prefix) + "accepted_raw_f64.bin").reshape((469, 100), order="F")
raw_tangent = np.fromfile(str(prefix) + "packed_tangent_f64.bin").reshape((469, 100), order="F")

def occupied_frame(values):
    """Mirror normalization, inactive projection, and active Cholesky gauge."""
    b = values / d[:, None]
    b /= np.sqrt(np.einsum("pi,pq,qi->i", b, s, b))[None, :]
    core = b[:, :89]
    rcore = np.linalg.cholesky(core.T @ s @ core).T
    core = solve_triangular(rcore.T, core.T, lower=True).T
    active = b[:, 89:] - core @ (core.T @ s @ b[:, 89:])
    ractive = np.linalg.cholesky(active.T @ s @ active).T
    active = solve_triangular(ractive.T, active.T, lower=True).T
    return np.column_stack((core, active))

occupied = occupied_frame(raw)
chol_s = np.linalg.cholesky(s)
white = chol_s.T @ occupied
virtual = solve_triangular(chol_s.T, null_space(white.T), lower=False)
mo = np.column_stack((occupied, virtual))
gram_error = float(np.max(np.abs(mo.T @ s @ mo - np.eye(469))))
assert gram_error < 1e-8
print("CROSS_FRAME", json.dumps(dict(ao_overlap_error=ao_error,
      gram_error=gram_error, available_cpus=len(os.sched_getaffinity(0)),
      pyscf_threads=lib.num_threads())), flush=True)

def masks(side):
    return [sum(1 << int(i) for i in line.split()) for line in
            Path(str(prefix) + side + "_strings.txt").read_text().splitlines()]
alpha, beta = masks("alpha"), masks("beta")
index_a, index_b = {m: i for i, m in enumerate(alpha)}, {m: i for i, m in enumerate(beta)}
ca = list(map(int, fci.cistring.make_strings(range(11), 5)))
cb = list(map(int, fci.cistring.make_strings(range(11), 3)))
ci = np.fromfile(str(prefix) + "carrier_f64.bin").reshape((462, 165), order="F")
ci = ci[np.ix_([index_a[m] for m in ca], [index_b[m] for m in cb])]
ci /= np.linalg.norm(ci)
mc = mcscf.CASSCF(mf, 11, (5, 3))
mc.max_memory = 120000
mc.conv_tol, mc.conv_tol_grad = 1e-7, 1e-3
mc.max_cycle_macro = 2
mc.fcisolver.conv_tol, mc.fcisolver.max_cycle = 1e-10, 200
mc.fcisolver.spin = 2
mc.fix_spin_(ss=2)
mc.chkfile = "cross-restart.chk"
mc.verbose = 5

# Compare physical directional gradients, not another fixed-point energy.
start = time.monotonic()
rdms = fci.direct_spin1.make_rdm12(ci, 11, (5, 3))
print("CROSS_RDM_READY", time.monotonic() - start, flush=True)
g = mc.get_grad(mo, rdms)
h = 1e-5
dm = (occupied_frame(raw + h * raw_tangent) -
      occupied_frame(raw - h * raw_tangent)) / (2 * h)
kappa = np.zeros((469, 469))
kappa[:, :100] = mo.T @ s @ dm
kappa[:100, 100:] = -kappa[100:, :100].T
print("CROSS_DIRECTION", json.dumps(dict(
    ao_overlap_error=ao_error, gram_error=gram_error,
    orbital_gradient_norm=float(np.linalg.norm(g)),
    # mc1step stores half the energy gradient; newton_casscf.gen_g_hop
    # explicitly returns g_orb * 2. Compare the energy derivative, not the
    # two optimizers' differently scaled internal covectors.
    mc1step_internal_contraction=float(g @ mc.pack_uniq_var(kappa)),
    directional_derivative=float(2 * g @ mc.pack_uniq_var(kappa)),
    occupied_tangent_orthogonality_error=float(np.linalg.norm(
        kappa[:100, :100] + kappa[:100, :100].T)),
    spin_square=float(fci.spin_op.spin_square(ci, 11, (5, 3))[0]))), flush=True)

if mode == "newton":
    mc = mc.newton()
elif mode != "default-ah":
    raise ValueError(mode)
mc.max_cycle_macro = 2
mc.kernel(mo, ci)
report = dict(mode=mode, converged=bool(mc.converged), energy=float(mc.e_tot),
              elapsed_seconds=time.monotonic() - start,
              spin_square=float(fci.spin_op.spin_square(mc.ci, 11, (5, 3))[0]))
Path("result.json").write_text(json.dumps(report, indent=2) + "\n")
print("CROSS_RESTART", json.dumps(report), flush=True)
