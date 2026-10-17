"""PySCF CASSCF references using the one shared AO/MO fixture."""

import argparse
import json
import os
from pathlib import Path
import time

if not os.environ.get("SLURM_JOB_ID"):
    raise SystemExit("molecular calculations require a Slurm compute allocation")

import numpy as np
import pyscf
from pyscf import fci, gto, lib, mcscf, scf
from pyscf.gto.basis import parse_gaussian
from pyscf.tools import molden

from fixture_format import read, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixture", type=Path)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("source_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--solver", choices=("newton", "ah"), default="newton",
                        help="coupled Newton CASSCF or standard orbital-only AH")
    parser.add_argument("--expected-case", choices=("F2", "FeCl2"),
                        help="reject a fixture for the wrong scheduled case")
    args = parser.parse_args()
    data = read(args.fixture, args.metadata)
    meta = data["metadata"]
    if args.expected_case and meta["case"] != args.expected_case:
        raise ValueError("fixture case differs from scheduled case")
    source = args.source_root.resolve()
    if (sha256(source / meta["basis"]) != meta["basis_sha256"] or
            sha256(source / meta["source_xmi"]) != meta["source_xmi_sha256"]):
        raise ValueError("source basis or molecular input changed after export")
    lib.num_threads(int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    basis_path = source / meta["basis"]
    atoms_bohr = [(name, [coord / 0.529177249 for coord in xyz])
                  for name, xyz in meta["atoms_angstrom"]]
    layout = meta["seed"].get("molden_ao_layout")
    if layout:
        layout = Path(layout)
        if sha256(layout) != meta["seed"]["molden_ao_layout_sha256"]:
            raise ValueError("qualified Molden AO layout changed after export")
        mol = molden.load(str(layout))[0]
        mol.charge = meta["charge"]
        mol.spin = data["spin"]
        mol.symmetry = False
        mol.set_geom_(np.asarray([xyz for _, xyz in atoms_bohr]), unit="Bohr",
                      symmetry=False, inplace=True)
        mol.max_memory, mol.verbose = 16000, 4
    else:
        basis = {name: parse_gaussian.load(str(basis_path), name, optimize=False)
                 for name in {atom[0] for atom in meta["atoms_angstrom"]}}
        mol = gto.M(atom=atoms_bohr, unit="Bohr", basis=basis, cart=True,
                    charge=meta["charge"], spin=data["spin"], symmetry=False,
                    max_memory=16000, verbose=4)
    if mol.nao_nr() != data["nbf"] or mol.has_ecp() or mol.nelectron != (
            2 * data["ncore"] + data["nelecas"]):
        raise ValueError("molecular identity differs from saved fixture")
    if abs(mol.energy_nuc() - data["nuclear"]) > 1e-9:
        raise ValueError("nuclear repulsion differs from saved fixture")
    mf = scf.RHF(mol) if data["spin"] == 0 else scf.ROHF(mol)
    # The saved fixture is the numerical authority for both solvers. Avoid
    # repeating SCF or regenerating ERIs and accidentally changing the input.
    mf.get_hcore = lambda *unused: data["hcore"]
    mf.get_ovlp = lambda *unused: data["overlap"]
    mf._eri = data["eri"]
    mf.mo_coeff = np.array(data["initial_mo"], order="F", copy=True)
    mf.mo_occ = np.asarray(meta["mo_occ"])
    mf.mo_energy = (np.asarray(meta["mo_energy"]) if meta["mo_energy"] is not None
                    else np.zeros(data["nbf"]))
    mf.e_tot = float(mf.energy_tot(dm=mf.make_rdm1(mf.mo_coeff, mf.mo_occ)))
    if layout:
        if abs(mf.e_tot - meta["seed"]["seed_determinant_energy"]) > 1e-8:
            raise ValueError("FeCl2 seed determinant energy differs from export")
    elif abs(mf.e_tot - meta["seed"]["hf_energy"]) > 1e-8:
        raise ValueError("F2 HF determinant energy differs from export")
    mf.converged = bool(meta["seed"]["hf_converged"])
    nalpha = (data["nelecas"] + data["spin"]) // 2
    nbeta = (data["nelecas"] - data["spin"]) // 2
    mc = mcscf.CASSCF(mf, data["ncas"], (nalpha, nbeta),
                      ncore=data["ncore"])
    mc.fcisolver.spin = data["spin"]
    mc.fix_spin_(ss=data["spin"] * (data["spin"] + 2) / 4, shift=0.5)
    # The Newton variant couples CI and orbital response. The ordinary
    # CASSCF object retains PySCF's standard orbital-only AH update.
    if args.solver == "newton":
        mc = mc.newton()
    mc.conv_tol = 1e-7
    mc.conv_tol_grad = 1e-3
    mc.max_cycle_macro = 100
    mc.max_memory = 16000
    mc.fcisolver.conv_tol = 1e-10
    mc.fcisolver.max_cycle = 200
    mc.verbose = 5
    trajectory = []

    def callback(env):
        entry = {"callback_index": len(trajectory) + 1,
                 "wall_seconds": time.monotonic() - start}
        for key in ("imacro", "e_tot", "e_cas", "de", "norm_gall",
                    "norm_gorb", "norm_gci", "tot_hop", "tot_kf"):
            value = env.get(key)
            if isinstance(value, (int, float, np.integer, np.floating)):
                entry[key] = float(value)
        stat = env.get("stat")
        if stat is not None:
            for key in ("imic", "tot_hop", "tot_kf"):
                value = getattr(stat, key, None)
                if isinstance(value, (int, float, np.integer, np.floating)):
                    entry["step_" + key] = float(value)
        trajectory.append(entry)
        print("CAS_MACRO " + json.dumps(entry, sort_keys=True), flush=True)
    mc.callback = callback

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    mc.kernel(mf.mo_coeff)
    last_hop = (trajectory[-1].get("tot_hop") if trajectory else None)
    result = dict(case=meta["case"],
                  solver=("PySCF newton_casscf coupled second order"
                          if args.solver == "newton" else
                          "PySCF standard CASSCF orbital-only AH"),
                  solver_mode=args.solver,
                  pyscf_version=pyscf.__version__,
                  fixture_sha256=meta["fixture_sha256"],
                  initial_casci_energy=meta["initial_casci_energy"],
                  converged=bool(mc.converged), energy=float(mc.e_tot),
                  elapsed_seconds=time.monotonic() - start,
                  macro_callback_count=len(trajectory), trajectory=trajectory,
                  last_reported_imacro=(int(trajectory[-1]["imacro"])
                                        if trajectory and "imacro" in trajectory[-1]
                                        else None),
                  hvp_actions=(int(last_hop) if last_hop is not None else None),
                  hvp_count_method=(
                      ("newton_casscf coupled kernel tot_hop" if args.solver == "newton"
                       else "standard CASSCF orbital-only AH kernel tot_hop")
                      if last_hop is not None else None),
                  final_spin_squared=float(fci.spin_op.spin_square(
                      mc.ci, mc.ncas, mc.nelecas)[0]),
                  native_threads=int(lib.num_threads()))
    (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print("CAS_RESULT " + json.dumps({k: v for k, v in result.items()
                                       if k != "trajectory"}, sort_keys=True), flush=True)
    if not mc.converged:
        raise RuntimeError("PySCF " + args.solver + " CASSCF did not converge")


if __name__ == "__main__":
    main()
