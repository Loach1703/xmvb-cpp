"""Compare saved orbital paths geometrically; no CI/energy reevaluations."""
import json
import os
from pathlib import Path
import re
import sys

assert os.environ.get("SLURM_JOB_ID"), "Use a Slurm compute node"
import numpy as np
import h5py
from scipy.linalg import solve_triangular
from pyscf import scf

root, run = Path(sys.argv[1]), Path(sys.argv[2])
mol, _ = scf.chkfile.load_scf(str(root / "porphine-reference.chk"))
deck = (root / "inputs/singlet/input.xmi").read_text()
geometry = []
for line in re.split(r"\$GEO\b", deck, flags=re.I)[1].splitlines():
    if line.strip().upper() == "$END":
        break
    if line.strip():
        fields = line.split()
        geometry.append((fields[0], tuple(float(x) / 0.529177249 for x in fields[1:])))
mol.spin = 0
mol.set_geom_(geometry, unit="Bohr", inplace=True)
mol.build(False, False)
# Only the one-electron overlap is evaluated to define a common geometry.
# No Hamiltonian, CI, gradient or Hessian evaluation is performed here.
s = mol.intor_symmetric("int1e_ovlp")
d = np.sqrt(s.diagonal())
initial = np.load(run / "pyscf-4/initial-mo.npy")[:, :100]


def frame(raw):
    b = raw / d[:, None]
    b /= np.sqrt(np.einsum("pi,pq,qi->i", b, s, b))[None, :]
    core = b[:, :89]
    upper = np.linalg.cholesky(core.T @ s @ core).T
    core = solve_triangular(upper.T, core.T, lower=True).T
    active = b[:, 89:] - core @ (core.T @ s @ b[:, 89:])
    upper = np.linalg.cholesky(active.T @ s @ active).T
    active = solve_triangular(upper.T, active.T, lower=True).T
    return np.column_stack((core, active))


def compare(left, right):
    result = {}
    for name, indices in (("core", slice(0, 89)), ("active", slice(89, 100)),
                          ("core_plus_active", slice(0, 100))):
        a, b = left[:, indices], right[:, indices]
        sigma = np.linalg.svd(a.T @ s @ b, compute_uv=False)
        sigma = np.clip(sigma, 0.0, 1.0)
        result[name] = dict(max_angle_degrees=float(np.degrees(np.arccos(sigma.min()))),
                            projector_distance=float(np.sqrt(max(0.0, 2*a.shape[1]-2*(sigma@sigma)))))
    result["left_core_right_active_weight"] = float(
        np.linalg.norm(left[:, :89].T @ s @ right[:, 89:])**2)
    return result


ao_locations = mol.ao_loc_nr()
fe_d_indices = [index for shell in range(mol.nbas)
                if mol.atom_pure_symbol(mol.bas_atom(shell)) == "Fe" and mol.bas_angular(shell) == 2
                for index in range(ao_locations[shell], ao_locations[shell+1])]
fe_gram = s[np.ix_(fe_d_indices, fe_d_indices)]


def character(mo):
    active = mo[:, 89:100]
    overlaps = s[fe_d_indices, :] @ active
    return float(np.trace(overlaps.T @ np.linalg.solve(fe_gram, overlaps)))


frames = {"initial": initial}
for count in (1, 4):
    raw = np.fromfile(run / f"neo-{count}/final-raw.bin").reshape((469, 100), order="F")
    frames[f"neo-{count}"] = frame(raw)
for path in sorted((run / "pyscf-4").glob("macro-*.npz")):
    with np.load(path) as saved:
        frames[f"pyscf-{int(saved['macro'])}"] = saved["mo"][:, :100]

# Read saved endpoints only. This is geometric path analysis, not another
# fixed-point Hamiltonian/gradient comparison or an optimization restart.
high_raw = np.fromfile(root / "porphine-terminal-1253506/cpp/audit_orbital_gradient_accepted_raw_f64.bin")
historical_high = frame(high_raw.reshape((469, 100), order="F"))
with h5py.File(root / "porphine-rotation-1248772/singlet/pyscf-ah/cross-restart.chk", "r") as file:
    reference_low = file["mcscf/mo_coeff"][:, :100]

report = dict(
    gram_errors={name: float(np.max(np.abs(mo.T@s@mo-np.eye(100)))) for name, mo in frames.items()},
    fe_cartesian_d_projector_weight={name: character(mo) for name, mo in frames.items()},
    from_initial={name: compare(initial, mo) for name, mo in frames.items() if name != "initial"},
    between_methods={str(count): compare(frames[f"neo-{count}"], frames[f"pyscf-{count}"])
                     for count in (1, 4)},
    to_saved_endpoints={name: dict(historical_high=compare(mo, historical_high),
                                  pyscf_low=compare(mo, reference_low))
                        for name, mo in frames.items()},
    note="Fe d projector includes all Cartesian Fe d AO functions, not an atomic charge or a pure Fe 3d/4d occupation.")


def printed_natural_density(path):
    """Use existing final NO output; record its six-decimal rounding error."""
    text = path.read_text().split("COMPUTED NATURAL ORBITALS", 1)[1]
    lines = text.splitlines()
    coefficients = np.zeros((469, 100))
    occupations = np.zeros(100)
    seen = np.zeros((469, 100), dtype=bool)
    row = 1
    while row < len(lines):
        fields = lines[row][21:].split()
        if not lines[row][:21].strip() and fields and all(re.fullmatch(r"\d+", x) for x in fields):
            columns = np.array([int(x)-1 for x in fields])
            row += 1
            occupations[columns] = list(map(float, lines[row][21:].split()))
            for _ in range(469):
                row += 1
                basis = int(lines[row][:21].split()[0])-1
                coefficients[basis, columns] = list(map(float, lines[row][21:].split()))
                seen[basis, columns] = True
        row += 1
        if seen.all():
            break
    assert seen.all(), path
    coefficients /= d[:, None]
    density = (coefficients * occupations) @ coefficients.T
    gram_error = float(np.max(np.abs(coefficients.T@s@coefficients-np.eye(100))))
    return density, occupations, gram_error


chol = np.linalg.cholesky(s)
for count in (1, 4):
    density, occupations, gram_error = printed_natural_density(run / f"neo-{count}/output.out")
    with np.load(run / f"pyscf-4/macro-{count:03d}.npz") as saved:
        py_mo = saved["mo"]
        py_density = 2*py_mo[:, :89] @ py_mo[:, :89].T
        py_density += py_mo[:, 89:100] @ saved["active_density"] @ py_mo[:, 89:100].T
    left = frames[f"neo-{count}"][:, 89:100]
    right = frames[f"pyscf-{count}"][:, 89:100]
    u, sigma, vt = np.linalg.svd(left.T@s@right)
    left_modes, right_modes = left@u, right@vt.T
    neo_mode_occupations = np.diag(left_modes.T@s@density@s@left_modes)
    py_mode_occupations = np.diag(right_modes.T@s@py_density@s@right_modes)
    def fe_weight(modes):
        overlaps = s[fe_d_indices, :] @ modes
        return np.diag(overlaps.T @ np.linalg.solve(fe_gram, overlaps))
    report.setdefault("active_principal_modes", {})[str(count)] = dict(
        angles_degrees=np.degrees(np.arccos(np.clip(sigma, 0, 1))).tolist(),
        neo_occupations=neo_mode_occupations.tolist(), pyscf_occupations=py_mode_occupations.tolist(),
        neo_fe_d_weights=fe_weight(left_modes).tolist(), pyscf_fe_d_weights=fe_weight(right_modes).tolist())
    report.setdefault("physical_density_comparison", {})[str(count)] = dict(
        frobenius_difference=float(np.linalg.norm(chol.T @ (density-py_density) @ chol)),
        neo_electron_trace=float(np.trace(density@s)),
        neo_printed_no_gram_error=gram_error,
        neo_printed_active_occupations=occupations[89:].tolist(),
        note="NEO natural-orbital coefficients and occupations are rounded to six decimals; this is not a derivative accuracy test.")
for count in range(1, 5):
    with np.load(run / f"pyscf-4/macro-{count:03d}.npz") as saved:
        report.setdefault("pyscf_active_occupations", {})[str(count)] = np.linalg.eigvalsh(saved["active_density"])[::-1].tolist()
if len(sys.argv) == 4:
    control = Path(sys.argv[3])
    with np.load(control / "macro-001.npz") as saved:
        control_mo = saved["mo"][:, :100]
        report["coupled_newton_control"] = dict(
            energy=float(saved["energy"]),
            from_initial=compare(initial, control_mo),
            vs_neo_first=compare(frames["neo-1"], control_mo),
            vs_default_ah_first=compare(frames["pyscf-1"], control_mo),
            fe_d_projector_weight=character(control_mo),
            active_occupations=np.linalg.eigvalsh(saved["active_density"])[::-1].tolist())
(run / "path-comparison.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2), flush=True)
