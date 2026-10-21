"""Bounded compute-node diagnostic for the archived FeCl2 Molden basis."""

import os
import sys

if not os.environ.get("SLURM_JOB_ID"):
    raise SystemExit("run under a Slurm compute allocation")

import numpy as np
from pyscf import gto
from pyscf.gto.basis import parse_gaussian
from pyscf.tools import molden


layout, basis_path = sys.argv[1:3]
atoms_bohr = [("Fe", (0., 0., 0.)),
              ("Cl", (0., 0., 2.136319 / 0.529177249)),
              ("Cl", (0., 0., -2.136319 / 0.529177249))]
expected_xyz = np.asarray([xyz for _, xyz in atoms_bohr])
mol = molden.load(layout)[0]
print("loaded", "nao", mol.nao_nr(), "cart", mol.cart,
      "natm", mol.natm, "symbols", [mol.atom_symbol(i) for i in range(mol.natm)],
      "units", mol.unit, flush=True)
print("loaded_xyz", mol.atom_coords().tolist(), flush=True)
mol.charge = 0
mol.spin = 4
mol.symmetry = False
mol.set_geom_(expected_xyz, unit="Bohr", symmetry=False, inplace=True)
print("set_geom", "nao", mol.nao_nr(), "cart", mol.cart,
      "natm", mol.natm, "symbols", [mol.atom_symbol(i) for i in range(mol.natm)],
      "xyz_max_diff", np.max(np.abs(mol.atom_coords() - expected_xyz)),
      flush=True)
mol.cart = True
mol.build(False, False)
print("forced_cart", "nao", mol.nao_nr(), "cart", mol.cart,
      "natm", mol.natm, "xyz_max_diff",
      np.max(np.abs(mol.atom_coords() - expected_xyz)), flush=True)
basis = {name: parse_gaussian.load(basis_path, name, optimize=False)
         for name in ("Fe", "Cl")}
named = gto.M(atom=atoms_bohr, unit="Bohr", basis=basis, cart=True,
              charge=0, spin=4, symmetry=False, verbose=0)
print("named", "nao", named.nao_nr(), "cart", named.cart,
      "labels_identical", mol.ao_labels() == named.ao_labels(), flush=True)
if mol.nao_nr() == named.nao_nr():
    print("overlap_max_diff", np.max(np.abs(
        mol.intor_symmetric("int1e_ovlp") - named.intor_symmetric("int1e_ovlp"))),
        flush=True)
