"""Geometry of saved porphine paths: no Hamiltonian, CI, energy or gradient calls."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re

assert os.environ.get("SLURM_JOB_ID"), "Use a Slurm compute node"
import h5py
import numpy as np
from scipy.linalg import expm, solve_triangular
from scipy.optimize import brentq
from pyscf import scf


def normalize(b, s):
    return b / np.sqrt(np.sum(b * (s @ b), axis=0))


def orthogonalize(b, s):
    r = np.linalg.cholesky(b.T @ s @ b).T
    return solve_triangular(r.T, b.T, lower=True).T, r


def geometry(b, s, n_core=89):
    b = normalize(b, s)
    qi, ri = orthogonalize(b[:, :n_core], s)
    a = b[:, n_core:]
    k = qi.T @ s @ a
    t = a - qi @ k
    qa, ra = orthogonalize(t, s)
    return dict(raw=b, core=qi, active=qa, occupied=np.column_stack((qi, qa)),
                ri=ri, ra=ra, k=k, active_gram=t.T @ s @ t)


def compare(left, right, s):
    result = {}
    for label in ("core", "active", "occupied"):
        sigma = np.linalg.svd(left[label].T @ s @ right[label], compute_uv=False)
        sigma = np.clip(sigma, 0, 1)
        result[label] = dict(angles_degrees=np.degrees(np.arccos(sigma)).tolist(),
                            max_angle_degrees=float(np.degrees(np.arccos(sigma.min()))),
                            projector_distance=float(np.sqrt(2 * np.sum(1 - sigma**2))))
    result["old_core_new_active_weight"] = float(np.linalg.norm(
        left["core"].T @ s @ right["active"]) ** 2)
    return result


def aligned_frame(old, new, s):
    u, _, vt = np.linalg.svd(old.T @ s @ new)
    return new @ vt.T @ u.T


def frame_rotation(old, new, s):
    """Oriented frame overlap is distinct from unsigned subspace angles."""
    u, _, vt = np.linalg.svd(old.T @ s @ new)
    polar = u @ vt
    return dict(determinant=float(np.linalg.det(polar)),
                distance_from_identity=float(np.linalg.norm(polar-np.eye(old.shape[1]))),
                max_angle_degrees=float(np.degrees(np.max(np.abs(
                    np.angle(np.linalg.eigvals(polar)))))))


def representative_decomposition(left, right, s):
    """Decompose normalized raw movement; do not label every core part gauge."""
    delta = right["raw"][:, 89:] - left["raw"][:, 89:]
    core = left["core"] @ (left["core"].T @ s @ delta)
    active = left["active"] @ (left["active"].T @ s @ delta)
    external = delta - core - active
    squared = lambda x: float(np.sum(x * (s @ x)))
    total = squared(delta)
    parts = dict(core_reference=squared(core), active_reference=squared(active),
                 external_reference=squared(external))
    # A geometry-only alternate representative of exactly the same new spans.
    # This is not a wavefunction/CI transformation for arbitrary selected VB.
    qi = aligned_frame(left["core"], right["core"], s)
    qa = aligned_frame(left["active"], right["active"], s)
    closest = np.column_stack((qi @ left["ri"], qi @ left["k"] + qa @ left["ra"]))
    cosine = np.sum(left["raw"] * (s @ right["raw"]), axis=0)
    cosine_aligned = np.sum(left["raw"] * (s @ closest), axis=0)
    return dict(raw_active_delta_squared=total, components_squared=parts,
                components_fraction={k: v / total if total else 0.0 for k, v in parts.items()},
                orthogonal_sum_error=abs(total - sum(parts.values())),
                actual_raw_active_angles_degrees=np.degrees(np.arccos(np.clip(cosine[89:], -1, 1))).tolist(),
                same_span_fixed_gram_angles_degrees=np.degrees(np.arccos(np.clip(cosine_aligned[89:], -1, 1))).tolist(),
                fixed_gram_error=float(np.max(np.abs(closest.T @ s @ closest -
                                                     left["raw"].T @ s @ left["raw"]))),
                same_span_check=compare(right, geometry(closest, s), s))


def read_guess(deck):
    block = re.split(r"\$GUS\b", deck.read_text(), flags=re.I)[1]
    raw = np.zeros((469, 100))
    orbital = None
    for line in block.splitlines():
        if line.strip().upper() == "$END":
            break
        match = re.search(r"#\s*ORBITAL\s+(\d+)\s+NAO", line, re.I)
        if match:
            orbital = int(match[1]) - 1
        elif orbital is not None:
            fields = line.split()
            for value, index in zip(fields[::2], fields[1::2]):
                raw[int(index)-1, orbital] = float(value.replace("D", "E"))
    return raw


def read_molden(path, mol):
    lines = path.read_text().split("[MO]", 1)[1].splitlines()
    columns = []
    current = None
    for line in lines:
        if line.strip().startswith("Sym="):
            if current is not None:
                columns.append(current)
            current = np.zeros(469)
        elif current is not None:
            fields = line.split()
            if len(fields) == 2 and fields[0].isdigit():
                current[int(fields[0])-1] = float(fields[1])
    if current is not None:
        columns.append(current)
    values = np.column_stack(columns)
    locations = mol.ao_loc_nr()
    permutation = []
    for shell in range(mol.nbas):
        count = locations[shell+1] - locations[shell]
        ell = mol.bas_angular(shell)
        order = {2: [0, 3, 5, 1, 2, 4], 3: [0, 6, 9, 3, 1, 2, 5, 8, 7, 4]}.get(ell, list(range(count)))
        assert len(order) == count and ell <= 3
        permutation.extend(locations[shell] + np.array(order))
    raw = np.empty_like(values)
    raw[np.array(permutation)] = values
    return raw[:, :100]


def tangent_diagnostics(origin, trial, s, metric_record):
    """Differentiate span projectors with the exact pre-canonicalization raw step."""
    i, a = origin[:, :89], origin[:, 89:]
    di, da = (trial-origin)[:, :89], (trial-origin)[:, 89:]
    gi = i.T @ s @ i
    qi, ri = orthogonalize(i, s)
    t = a - qi @ (qi.T @ s @ a)
    qa, ra = orthogonalize(t, s)
    h = di - qi @ (qi.T @ s @ di)
    dt = da - qi @ (qi.T @ s @ da)
    dt -= h @ np.linalg.solve(gi, i.T @ s @ a)
    dt -= i @ np.linalg.solve(gi, h.T @ s @ a)
    z = solve_triangular(ri.T, h.T, lower=True).T
    y = dt - qa @ (qa.T @ s @ dt)
    y = solve_triangular(ra.T, y.T, lower=True).T
    norm2 = lambda x: float(np.sum(x * (s @ x)))
    independent = np.sqrt(norm2(z) + norm2(y))
    old, new = geometry(origin, s), geometry(trial, s)
    upper = np.linalg.cholesky(s).T
    result = dict(independent_physical_tangent_norm=float(independent),
                  logged_physical_tangent_norm=metric_record,
                  geometry=compare(old, new, s),
                  raw=representative_decomposition(old, new, s), blocks={})
    if metric_record is not None:
        result["metric_norm_error"] = abs(float(independent) - metric_record)
    # Cholesky orthogonal frames can rotate inside the active span even when
    # the raw tangent is orthogonal to gauge directions in Euclidean AO units.
    dt_frame = solve_triangular(ra.T, dt.T, lower=True).T
    parallel = qa.T @ s @ dt_frame
    chol_connection = np.triu(parallel+parallel.T)
    chol_connection[np.diag_indices_from(chol_connection)] *= 0.5
    skew_connection = parallel-chol_connection
    singular = np.linalg.svd(skew_connection, compute_uv=False)
    result["active_frame_tangent"] = dict(
        parallel_frobenius=float(np.linalg.norm(parallel)),
        cholesky_connection_frobenius=float(np.linalg.norm(chol_connection)),
        rotation_frobenius=float(np.linalg.norm(skew_connection)),
        rotation_max_radians=float(singular.max()),
        rotation_skew_error=float(np.linalg.norm(skew_connection+skew_connection.T)),
        physical_external_frobenius=float(np.sqrt(norm2(y))))
    new_t = trial[:, 89:] - new["core"] @ (new["core"].T @ s @ trial[:, 89:])
    w = solve_triangular(ra.T, new_t.T, lower=True).T
    finite_parallel = qa.T @ s @ w
    finite_external = w-qa @ finite_parallel
    sv = np.linalg.svd(finite_parallel, compute_uv=False)
    graph = np.linalg.solve(finite_parallel.T, finite_external.T).T
    result["active_finite_graph"] = dict(
        parallel_singular_values=sv.tolist(),
        parallel_condition=float(sv.max()/sv.min()),
        gram_condition=float(np.linalg.cond(new_t.T @ s @ new_t)),
        external_frobenius=float(np.sqrt(norm2(finite_external))),
        graph_frobenius=float(np.sqrt(norm2(graph))))
    result["active_finite_frame_rotation"] = frame_rotation(qa, new["active"], s)
    if result["active_finite_frame_rotation"]["determinant"] < 0:
        def determinant(scale):
            point = geometry(origin+scale*(trial-origin), s)
            return np.linalg.det(qa.T @ s @ point["active"])
        crossing = brentq(determinant, 0, 1, xtol=1e-12)
        point = geometry(origin+crossing*(trial-origin), s)
        result["oriented_frame_crossing"] = dict(
            fraction=float(crossing), geometry=compare(old, point, s),
            active_gram_condition=float(np.linalg.cond(point["active_gram"])))
    for label, derivative in (("core", z), ("active", y)):
        q0, q1, dw = upper @ old[label], upper @ new[label], upper @ derivative
        linear = dw @ q0.T + q0 @ dw.T
        actual = q1 @ q1.T - q0 @ q0.T
        result["blocks"][label] = dict(
            linear_projector_norm=float(np.linalg.norm(linear)),
            actual_projector_norm=float(np.linalg.norm(actual)),
            nonlinear_remainder_norm=float(np.linalg.norm(actual-linear)))
    # Geometry-only isometry with the SAME first physical tangent, not a new
    # optimizer or a claim that arbitrary selected VB structures are closed.
    q = upper @ np.column_stack((qi, qa))
    d = upper @ np.column_stack((z, y))
    connection = q.T @ d
    omega = d @ q.T - q @ d.T - q @ connection @ q.T
    result["isometry_tangent"] = dict(
        connection_skew_error=float(np.linalg.norm(connection + connection.T)),
        generator_skew_error=float(np.linalg.norm(omega + omega.T)),
        reconstruction_error=float(np.linalg.norm(omega @ q-d)))
    q_exp = expm(omega) @ q
    back = solve_triangular(upper, q_exp, lower=False)
    iso = dict(core=back[:, :89], active=back[:, 89:], occupied=back)
    result["same_tangent_isometry_geometry"] = compare(old, iso, s)
    result["raw_vs_same_tangent_isometry"] = compare(new, iso, s)
    # Independent local spheres have the same normalized raw first tangent;
    # unlike the isometry they need not preserve any inter-orbital Gram.
    norms = np.sqrt(np.sum(origin*(s @ origin), axis=0))
    x = origin/norms
    d_raw = (trial-origin)/norms
    alpha = np.sum(x*(s @ d_raw), axis=0)
    tangent = d_raw-x*alpha
    theta = np.sqrt(np.maximum(0, np.sum(tangent*(s @ tangent), axis=0)))
    sphere = x*np.cos(theta)+tangent*np.sinc(theta/np.pi)
    sphere_frame = geometry(sphere, s)
    result["local_sphere_geometry"] = compare(old, sphere_frame, s)
    result["raw_vs_local_sphere"] = compare(new, sphere_frame, s)
    result["local_sphere_frame_rotation"] = frame_rotation(old["active"], sphere_frame["active"], s)
    for label in ("core", "active"):
        q0, q1 = upper @ old[label], upper @ iso[label]
        derivative = upper @ (z if label == "core" else y)
        linear = derivative @ q0.T + q0 @ derivative.T
        result["blocks"][label]["isometry_nonlinear_remainder_norm"] = float(
            np.linalg.norm(q1 @ q1.T-q0 @ q0.T-linear))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--snapshots", type=Path)
    parser.add_argument("--terminal", type=Path,
                        help="Additional native terminal Molden file to compare with the saved PySCF endpoint")
    args = parser.parse_args()
    root = args.root
    mol, _ = scf.chkfile.load_scf(str(root / "porphine-reference.chk"))
    deck = root / "inputs/singlet/input.xmi"
    coordinates = []
    for line in re.split(r"\$GEO\b", deck.read_text(), flags=re.I)[1].splitlines():
        if line.strip().upper() == "$END":
            break
        fields = line.split()
        if fields:
            coordinates.append((fields[0], tuple(float(x)/0.529177249 for x in fields[1:])))
    mol.spin = 0
    mol.set_geom_(coordinates, unit="Bohr", inplace=True)
    mol.build(False, False)
    # One-electron overlap only. No electronic Hamiltonian or many-body action.
    py_s = mol.intor_symmetric("int1e_ovlp")
    scale = np.sqrt(py_s.diagonal())
    s = py_s / scale[:, None] / scale[None, :]
    exported_s = np.fromfile(root / "porphine-coupled-audit-1247963/audit_orbital_gradient_ao_overlap_f64.bin")
    exported_s = exported_s.reshape((469, 469), order="F")
    assert np.max(np.abs(s-exported_s)) < 1e-8
    raw = dict(initial=read_guess(deck),
               neo1=read_molden(root / "transport-regression-1254179/singlet/input.molden", mol),
               neo2=read_molden(root / "step-geometry-1254271/singlet/input.molden", mol),
               neo_terminal=read_molden(root / "porphine-transport-1254189/singlet/input.molden", mol))
    if args.terminal:
        raw["neo_current_terminal"] = read_molden(args.terminal, mol)
    if args.snapshots:
        # Use the audited run's accepted frames, not stale pre-fix endpoints.
        raw["neo1"] = np.fromfile(args.snapshots / "macro-002-origin.bin").reshape((469, 100), order="F")
        raw["neo2"] = read_molden(args.snapshots.parent / "input.molden", mol)
    py_path = root / "porphine-path-1253884/pyscf-4"
    raw["pyscf_initial"] = np.load(py_path / "initial-mo.npy")[:, :100] * scale[:, None]
    for macro in range(1, 5):
        with np.load(py_path / f"macro-{macro:03d}.npz") as saved:
            raw[f"pyscf{macro}"] = saved["mo"][:, :100] * scale[:, None]
    with h5py.File(root / "porphine-rotation-1248772/singlet/pyscf-ah/cross-restart.chk") as checkpoint:
        raw["pyscf_low"] = checkpoint["mcscf/mo_coeff"][:, :100] * scale[:, None]
    exact = np.fromfile(root / "porphine-path-1253884/neo-1/final-raw.bin").reshape((469, 100), order="F")
    rounded = read_molden(root / "porphine-path-1253884/neo-1/input.molden", mol)
    frames = {name: geometry(b, s) for name, b in raw.items()}
    report = dict(provenance=dict(slurm_job_id=os.environ["SLURM_JOB_ID"],
                                 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                                 snapshots=str(args.snapshots) if args.snapshots else None,
                                 current_terminal=str(args.terminal) if args.terminal else None,
                                 current_terminal_sha256=hashlib.sha256(args.terminal.read_bytes()).hexdigest()
                                 if args.terminal else None,
                                 baseline_terminal_native_job="1254189"),
                  ao_overlap_error=float(np.max(np.abs(s-exported_s))),
                  molden_coefficient_rounding=float(np.max(np.abs(exact-rounded))),
                  molden_geometry_rounding=compare(geometry(exact, s), geometry(rounded, s), s),
                  initial_match=compare(frames["initial"], frames["pyscf_initial"], s),
                  frames={name: dict(active_gram_condition=float(np.linalg.cond(frame["active_gram"])))
                          for name, frame in frames.items()}, from_initial={}, between={}, consecutive={})
    for name, frame in frames.items():
        report["from_initial"][name] = dict(geometry=compare(frames["initial"], frame, s),
                                         raw=representative_decomposition(frames["initial"], frame, s))
    for left, right in (("neo1", "pyscf1"), ("neo2", "pyscf2"),
                        ("neo_terminal", "pyscf_low")):
        report["between"][left+"_"+right] = compare(frames[left], frames[right], s)
    if args.terminal:
        for right in ("pyscf_low", "neo_terminal"):
            report["between"]["neo_current_terminal_"+right] = compare(
                frames["neo_current_terminal"], frames[right], s)
    for left, right in (("neo1", "neo2"), ("pyscf1", "pyscf2"),
                        ("pyscf2", "pyscf3"), ("pyscf3", "pyscf4")):
        report["consecutive"][left+"_"+right] = dict(
            geometry=compare(frames[left], frames[right], s),
            raw=representative_decomposition(frames[left], frames[right], s),
            active_frame_rotation=frame_rotation(frames[left]["active"], frames[right]["active"], s))
    if args.snapshots:
        import csv
        native = args.snapshots.parent
        metric = {int(row["iteration"]): float(row["orbital_step_norm"])
                  for row in csv.DictReader((native / "neo.tsv").open(), delimiter="\t")}
        report["candidate_geometry"] = {}
        for path in sorted(args.snapshots.glob("macro-*.bin")):
            if path.name.endswith("-origin.bin"):
                continue
            macro = int(path.name.split("-")[1])
            origin = np.fromfile(args.snapshots / f"macro-{macro:03d}-origin.bin").reshape((469, 100), order="F")
            trial = np.fromfile(path).reshape((469, 100), order="F")
            result = tangent_diagnostics(origin, trial, s,
                                         metric.get(macro) if "macro_accepted" in path.name else None)
            if "macro_accepted" in path.name and macro == max(metric):
                canonical = geometry(read_molden(native / "input.molden", mol), s)
                result["canonicalization_geometry"] = compare(geometry(trial, s), canonical, s)
            report["candidate_geometry"][path.name] = result
            result["snapshot_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "representative-analysis.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
