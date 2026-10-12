"""Inspect saved basis/ECP metadata and measure spin of exported CI carriers.

No SCF/CASCI or energy comparison is performed.
"""
import hashlib
import json
import os
from pathlib import Path
import sys

assert os.environ.get("SLURM_JOB_ID"), "Use a Slurm compute node"
import numpy as np
from pyscf import fci, scf
from pyscf.gto.basis import parse_gaussian

root = Path(sys.argv[1])
mol, _ = scf.chkfile.load_scf(str(root / "porphine-reference.chk"))
basis_file = root / "inputs/triplet/cc-pvdz.gbs"
symbols = sorted(set(mol.atom_symbol(i) for i in range(mol.natm)))
comparison = {}
for symbol in symbols:
    expected = parse_gaussian.load(str(basis_file), symbol, optimize=False)
    comparison[symbol] = dict(
        identical=expected == mol._basis[symbol],
        shells=len(expected),
        primitives=sum(len(shell) - 1 for shell in expected))
report = dict(basis_sha256=hashlib.sha256(basis_file.read_bytes()).hexdigest(),
              per_element=comparison, ecp=mol._ecp, pseudo=mol.pseudo,
              nuclear_charges=mol.atom_charges().tolist(),
              electrons=mol.nelectron, charge=mol.charge, cartesian=mol.cart,
              n_bf=mol.nao_nr(), n_atoms=mol.natm, symmetry=mol.symmetry)
print("BASIS_AUDIT", json.dumps(report), flush=True)
assert all(item["identical"] for item in comparison.values())
assert not mol._ecp and not mol.pseudo and mol.nelectron == 186

if len(sys.argv) > 2:
    artifact = Path(sys.argv[2])
    prefix = artifact / "audit_orbital_gradient_"
    def strings(side):
        rows = Path(str(prefix) + side + "_strings.txt").read_text().splitlines()
        return [sum(1 << int(i) for i in line.split()) for line in rows]
    alpha, beta = strings("alpha"), strings("beta")
    na, nb = bin(alpha[0]).count("1"), bin(beta[0]).count("1")
    n_orb = 11
    canonical_a = list(map(int, fci.cistring.make_strings(range(n_orb), na)))
    canonical_b = list(map(int, fci.cistring.make_strings(range(n_orb), nb)))
    assert set(alpha) == set(canonical_a) and set(beta) == set(canonical_b)
    carrier = np.fromfile(str(prefix) + "carrier_f64.bin", dtype=np.float64)
    carrier = carrier.reshape((len(alpha), len(beta)), order="F")
    index_a, index_b = {m: i for i, m in enumerate(alpha)}, {m: i for i, m in enumerate(beta)}
    ci = carrier[np.ix_([index_a[m] for m in canonical_a],
                        [index_b[m] for m in canonical_b])]
    ci /= np.linalg.norm(ci)
    ss, multiplicity = fci.spin_op.spin_square(ci, n_orb, (na, nb))
    ss_image = fci.spin_op.contract_ss(ci, n_orb, (na, nb))
    target = ((na - nb) / 2) * ((na - nb) / 2 + 1)
    print("CARRIER_SPIN_AUDIT", json.dumps(dict(
        alpha_strings=len(alpha), beta_strings=len(beta),
        active_electrons=[na, nb], spin_square=float(ss),
        multiplicity=float(multiplicity), target=target,
        spin_eigenvector_residual=float(np.linalg.norm(ss_image - target * ci)))), flush=True)
