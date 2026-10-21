"""Create one exact AO and initial-MO fixture per case inside a Slurm job."""

import argparse
import json
import os
from pathlib import Path

if not os.environ.get("SLURM_JOB_ID"):
    raise SystemExit("molecular calculations require a Slurm compute allocation")

import numpy as np
import pyscf
from pyscf import fci, gto, lib, mcscf, scf, symm
from pyscf.gto.basis import parse_gaussian
from pyscf.tools import molden

from fecl2_seed import complete_orthonormal_mos, read_xmvb_orbitals
from fixture_format import HEADER, MAGIC, sha256


CASES = {
    "F2": dict(atoms=[("F", (0., 0., 0.)), ("F", (0., 0., 1.4))],
               basis="basis/cc-pvdz.gbs", xmi="testdata/vbscf/F2_OEO_PYSCF_LOW.xmi",
               charge=0, spin=0, nelecas=2, ncas=2),
    "FeCl2": dict(atoms=[("Fe", (0., 0., 0.)),
                          ("Cl", (0., 0., 2.136319)),
                          ("Cl", (0., 0., -2.136319))],
                  basis="basis/def2-tzvp.gbs", xmi="testdata/vbscf/FeCl2.xmi",
                  charge=0, spin=4, nelecas=10, ncas=8),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=CASES)
    parser.add_argument("source_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--fecl2-molden-layout", type=Path,
                        help="qualified XMVB input.molden AO shell layout; required for FeCl2")
    args = parser.parse_args()
    spec = CASES[args.case]
    source = args.source_root.resolve()
    basis_path = source / spec["basis"]
    xmi_path = source / spec["xmi"]
    lib.num_threads(int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    basis = {name: parse_gaussian.load(str(basis_path), name, optimize=False)
             for name in {atom[0] for atom in spec["atoms"]}}
    # Match xmvb's Angstrom-to-Bohr conversion used for Gaussian inputs.
    atoms_bohr = [(name, [coord / 0.529177249 for coord in xyz])
                  for name, xyz in spec["atoms"]]
    layout_path = None
    if args.case == "FeCl2":
        if args.fecl2_molden_layout is None:
            parser.error("FeCl2 requires --fecl2-molden-layout")
        layout_path = args.fecl2_molden_layout.resolve(strict=True)
        # Molden preserves the qualified XMVB AO shell order. Named PySCF
        # bases can group shells differently, invalidating $GUS coefficients.
        mol = molden.load(str(layout_path))[0]
        mol.charge = spec["charge"]
        mol.spin = spec["spin"]
        mol.symmetry = False
        mol.set_geom_(np.asarray([xyz for _, xyz in atoms_bohr]), unit="Bohr",
                      symmetry=False, inplace=True)
        mol.max_memory = 16000
        mol.verbose = 4
        # Molden.load retains indexed atom labels (Fe1, Cl2, Cl3); compare
        # physical nuclear charges instead of those display labels.
        if (mol.nao_nr() != 136 or not mol.cart or mol.natm != 3 or
                [mol.atom_charge(i) for i in range(3)] != [26, 17, 17] or
                np.max(np.abs(mol.atom_coords() - np.asarray([xyz for _, xyz in atoms_bohr]))) > 1e-8):
            raise RuntimeError("qualified FeCl2 Molden AO layout is not Cartesian 136-AO")
    else:
        mol = gto.M(atom=atoms_bohr, unit="Bohr", basis=basis, cart=True,
                    charge=spec["charge"], spin=spec["spin"], symmetry=False,
                    max_memory=16000, verbose=4)
    if mol.has_ecp():
        raise RuntimeError("ECP is outside this exact all-electron benchmark")
    ncore = (mol.nelectron - spec["nelecas"]) // 2
    if 2 * ncore + spec["nelecas"] != mol.nelectron:
        raise RuntimeError("electron count is inconsistent with inactive core")
    nbf = mol.nao_nr()
    mf = scf.RHF(mol) if args.case == "F2" else scf.ROHF(mol)
    s = np.ascontiguousarray(mf.get_ovlp(), dtype="<f8")
    h = np.ascontiguousarray(mf.get_hcore(), dtype="<f8")
    if args.case == "F2":
        mf.conv_tol = 1e-10
        mf.max_cycle = 200
        mf.kernel()
        if not mf.converged:
            raise RuntimeError("F2 RHF guess did not converge")
        sym_mol = gto.M(atom=atoms_bohr, unit="Bohr", basis=basis, cart=True,
                        charge=0, spin=0, symmetry="D2h", verbose=0)
        if (sym_mol.ao_labels() != mol.ao_labels() or
                np.max(np.abs(sym_mol.intor_symmetric("int1e_ovlp") - s)) > 1e-10):
            raise RuntimeError("F2 symmetry labeling changed the AO basis")
        labels = np.asarray(symm.label_orb_symm(
            sym_mol, sym_mol.irrep_name, sym_mol.symm_orb,
            mf.mo_coeff, s=s, check=False))
        occupied = np.flatnonzero(mf.mo_occ == 2)
        virtual = np.flatnonzero(mf.mo_occ == 0)
        ag = [i for i in occupied if labels[i] == "Ag"]
        b1u = [i for i in virtual if labels[i] == "B1u"]
        if not ag or not b1u:
            raise RuntimeError("F2 lacks occupied Ag or virtual B1u orbital")
        sigma_g = max(ag, key=lambda i: mf.mo_energy[i])
        sigma_u = min(b1u, key=lambda i: mf.mo_energy[i])
        def irrep_weight(irrep, column):
            block = sym_mol.symm_orb[sym_mol.irrep_name.index(irrep)]
            projected = block.T @ s @ mf.mo_coeff[:, column]
            metric = block.T @ s @ block
            return float(projected @ np.linalg.solve(metric, projected))
        weights = (irrep_weight("Ag", sigma_g), irrep_weight("B1u", sigma_u))
        if min(weights) < 0.99 or len(occupied) - 1 != ncore:
            raise RuntimeError("F2 sigma orbitals lack pure D2h character")
        order = np.r_[occupied[occupied != sigma_g], sigma_g, sigma_u,
                      virtual[virtual != sigma_u]]
        u = np.ascontiguousarray(mf.mo_coeff[:, order], dtype="<f8")
        mo_occ = np.asarray(mf.mo_occ[order], dtype=float)
        mo_energy = np.asarray(mf.mo_energy[order], dtype=float)
        seed = dict(kind="RHF F2 sigma_g/sigma_u_star selection",
                    occupied_ag_original_index=int(sigma_g),
                    virtual_b1u_original_index=int(sigma_u),
                    selected_irrep_weights=weights,
                    hf_energy=float(mf.e_tot), hf_converged=True)
    else:
        guess = read_xmvb_orbitals(xmi_path, nbf, ncore + spec["ncas"],
                                   np.sqrt(np.diag(s)))
        u = np.ascontiguousarray(complete_orthonormal_mos(s, guess, ncore),
                                 dtype="<f8")
        mo_occ = np.r_[np.full(ncore + 3, 2.), np.full(4, 1.),
                       np.zeros(nbf - ncore - 7)]
        mo_energy = None
        mf.mo_coeff = u
        mf.mo_occ = mo_occ
        seed = dict(kind="FeCl2 archived XMVB $GUS plus overlap-metric QR completion",
                    molden_ao_layout=str(layout_path),
                    molden_ao_layout_sha256=sha256(layout_path),
                    archive_seed_columns=ncore + spec["ncas"],
                    hf_energy=None, hf_converged=False)
    gram_error = float(np.max(np.abs(u.T @ s @ u - np.eye(nbf))))
    if u.shape != (nbf, nbf) or gram_error > 1e-8:
        raise RuntimeError("shared MO frame is not complete and S-orthonormal")
    nalpha = (spec["nelecas"] + spec["spin"]) // 2
    nbeta = (spec["nelecas"] - spec["spin"]) // 2
    if not (0 <= nbeta <= nalpha <= spec["ncas"]):
        raise RuntimeError("invalid CAS electrons or spin")
    # Compute and retain exactly one AO ERI stream; the initial PySCF CASCI
    # and both later CASSCF solvers consume this numerical source.
    eri = mol.intor("int2e", aosym="s8")
    npair = nbf * (nbf + 1) // 2
    neri = npair * (npair + 1) // 2
    if eri.shape != (neri,) or not np.isfinite(eri).all():
        raise RuntimeError("unexpected exact AO ERI shape or value")
    mf._eri = eri
    if args.case == "FeCl2":
        seed["seed_determinant_energy"] = float(
            mf.energy_tot(dm=mf.make_rdm1(u, mo_occ)))
    cas = mcscf.CASCI(mf, spec["ncas"], (nalpha, nbeta), ncore=ncore)
    cas.fcisolver.spin = spec["spin"]
    cas.fix_spin_(ss=spec["spin"] * (spec["spin"] + 2) / 4, shift=0.5)
    cas.fcisolver.conv_tol = 1e-11
    cas.kernel(u)
    if not cas.converged:
        raise RuntimeError("initial CASCI failed; no fixture exported")
    initial_casci = float(cas.e_tot)
    initial_spin_square = float(fci.spin_op.spin_square(
        cas.ci, spec["ncas"], cas.nelecas)[0])
    target_spin_square = spec["spin"] * (spec["spin"] + 2) / 4
    if abs(initial_spin_square - target_spin_square) > 1e-6:
        raise RuntimeError("initial CASCI root has the wrong total spin")
    # PySCF s8 is the lower triangle of the lower AO-pair matrix; this is
    # exactly xmvb's i*(i+1)/2+j pair-of-pairs order.
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    fixture = output / "ao-mo.bin"
    with fixture.open("wb") as stream:
        stream.write(HEADER.pack(MAGIC, nbf, ncore, spec["ncas"],
                                 spec["nelecas"], spec["spin"], mol.energy_nuc()))
        s.tofile(stream)
        h.tofile(stream)
        u.tofile(stream)
        np.asarray(eri, dtype="<f8").tofile(stream)
    metadata = dict(case=args.case, atoms_angstrom=spec["atoms"], basis=spec["basis"],
                    source_xmi=spec["xmi"], source_xmi_sha256=sha256(xmi_path),
                    basis_sha256=sha256(basis_path), fixture_sha256=sha256(fixture),
                    nbf=nbf, ncore=ncore, ncas=spec["ncas"],
                    nelecas=spec["nelecas"], spin=spec["spin"],
                    charge=spec["charge"], cartesian=True, symmetry=False,
                    integral_approximation="exact s8; no RI; no ECP",
                    ao_basis_contract=("qualified XMVB Molden shell order" if layout_path
                                       else "PySCF GBS shell order"),
                    seed=seed, mo_occ=mo_occ.tolist(),
                    mo_energy=(mo_energy.tolist() if mo_energy is not None else None),
                    initial_casci_energy=initial_casci,
                    initial_spin_square=initial_spin_square,
                    initial_gram_max_error=gram_error, pyscf_version=pyscf.__version__)
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
