"""Prepare and certify Molcas's AO conversion; run only on Slurm."""
import json
import os
from pathlib import Path
import sys

assert os.environ.get("SLURM_JOB_ID")
import h5py
import numpy as np
from pyscf import scf

root, stage = Path(sys.argv[1]), sys.argv[2]
mol, _ = scf.chkfile.load_scf(str(root / "porphine-reference.chk"))
assert mol.cart and not mol._ecp and mol.nelectron == 186
if stage == "prepare":
    names = {
        "H": "H.cc-pVDZ.Dunning.4s1p.2s1p.",
        "C": "C.cc-pVDZ.Dunning.9s4p1d.3s2p1d.",
        "N": "N.cc-pVDZ.Dunning.9s4p1d.3s2p1d.",
        "Fe": "Fe.cc-pVDZ.Balabanov.20s16p8d2f.6s5p3d1f.",
    }
    lines = ["&GATEWAY", "Title = Fe porphine triplet shared low-energy frame", "Group = Nosym"]
    for atom in range(mol.natm):
        symbol = mol.atom_pure_symbol(atom)
        xyz = " ".join(f"{x:.16e}" for x in mol.atom_coord(atom))
        lines += ["Basis set", names[symbol], f"{symbol}{atom+1} {xyz} /bohr",
                  "Cartesian all", "End of basis"]
    # Use conventional integrals, not Cholesky/RI or a relativistic/ECP Hamiltonian.
    gateway = "\n".join(lines) + "\n&SEWARD\n"
    Path("setup.input").write_text(gateway + "&GUESSORB\n")
    rasscf = """&RASSCF
Title = Fe porphine CAS(8,11) triplet from shared low-energy frame
FileOrb = $CurrDir/initial.RasOrb
Spin = 3
Charge = 0
Nactel = 8 0 0
Inactive = 89
Ras2 = 11
Thrs = 1.0e-7 1.0e-3 1.0e-3
Tight = 1.0e-10 1.0e-4
Iterations = 1000 200
Peri
QUNE
"""
    # Reuse the gateway/seward files already created by setup in the same scratch.
    Path("optimize.input").write_text(rasscf)
    sys.exit(0)
if stage != "convert":
    raise ValueError(stage)

with h5py.File("setup.guessorb.h5", "r") as f:
    assert int(f.attrs["NSYM"]) == 1
    ids = f["BASIS_FUNCTION_IDS"][:]
    n = len(ids)
    sm = f["AO_OVERLAP_MATRIX"][:].reshape(n, n)
    primitive_ids, primitives = f["PRIMITIVE_IDS"][:], f["PRIMITIVES"][:]
    centers, charges = f["CENTER_COORDINATES"][:], f["CENTER_CHARGES"][:]
assert n == 469 and ids.shape == (n, 4)
center_map = {}
for center, (xyz, charge) in enumerate(zip(centers, charges), 1):
    matches = [i for i in range(mol.natm) if mol.atom_charge(i) == charge
               and np.max(np.abs(mol.atom_coord(i) - xyz)) < 1e-9]
    assert len(matches) == 1, (center, matches)
    center_map[center] = matches[0]
assert sorted(center_map.values()) == list(range(mol.natm))

locations = mol.ao_loc_nr()
components, phases = {}, {}
for center in sorted(set(ids[:, 0])):
    for l in sorted(set(abs(int(row[2])) for row in ids if row[0] == center)):
        for radial in sorted(set(int(row[1]) for row in ids
                                 if row[0] == center and abs(row[2]) == l)):
            rows = primitives[np.all(primitive_ids == [center, l, radial], axis=1)]
            rows = rows[rows[:, 1] != 0]
            rows = rows[np.argsort(-rows[:, 0])]
            direction = rows[:, 1] / np.linalg.norm(rows[:, 1])
            matches = []
            for shell in range(mol.nbas):
                if mol.bas_atom(shell) != center_map[center] or mol.bas_angular(shell) != l:
                    continue
                for column, coeff in enumerate(mol.bas_ctr_coeff(shell).T):
                    mask = coeff != 0
                    exps = mol.bas_exp(shell)[mask]
                    if len(exps) != len(rows) or not np.allclose(exps, rows[:, 0], rtol=1e-11, atol=1e-12):
                        continue
                    candidate = coeff[mask] / np.linalg.norm(coeff[mask])
                    phase = 1 if candidate @ direction > 0 else -1
                    if np.max(np.abs(candidate - phase * direction)) < 1e-9:
                        matches.append((int(locations[shell]) + column * (l+1)*(l+2)//2, phase))
            assert len(matches) == 1, (center, l, radial, matches)
            offset, phase = matches[0]
            phases[(center, radial, l)] = phase
            component = 0
            for x in range(l, -1, -1):
                for y in range(l-x, -1, -1):
                    components[(center, radial, l, x, y, l-x-y)] = offset + component
                    component += 1
mapping, signs = [], []
for center, radial, angular, m in ids:
    if angular == 0:
        l, x, y, z = 0, 0, 0, 0
    elif angular == 1:
        l = 1
        x, y, z = {1: (1, 0, 0), -1: (0, 1, 0), 0: (0, 0, 1)}[int(m)]
    else:
        assert angular < 0, "Spherical shell detected"
        l = -int(angular)
        yz = int(np.floor((np.sqrt(8 * (int(m)+l)+1)-1)/2))
        z = int(m)+l-yz*(yz+1)//2
        x, y = l-yz, yz-z
    mapping.append(components[(center, radial, l, x, y, z)])
    signs.append(phases[(center, radial, l)])
assert sorted(mapping) == list(range(n))
sp = mol.intor_symmetric("int1e_ovlp")
dp, dm = np.sqrt(sp.diagonal()), np.sqrt(sm.diagonal())
signs = np.asarray(signs)
overlap_error = float(np.max(np.abs(sm / np.outer(dm, dm) -
    (sp / np.outer(dp, dp))[np.ix_(mapping, mapping)] * np.outer(signs, signs))))
assert overlap_error < 1e-8, overlap_error
mo = np.load("shared-mo.npy")
converted = mo[mapping, :] * (signs * dp[mapping] / dm)[:, None]
orth_error = float(np.max(np.abs(converted.T @ sm @ converted - np.eye(n))))
assert orth_error < 1e-8, orth_error
lines = ["#INPORB 2.2", "#INFO", "* Certified shared CPP/PySCF triplet frame",
         "       0       1       0", f"{n:8d}", f"{n:8d}", "#ORB"]
for orbital in range(n):
    lines.append(f"* ORBITAL    1 {orbital+1:5d}")
    for first in range(0, n, 5):
        lines.append("".join(f" {x:21.14E}" for x in converted[first:first+5, orbital]))
occupation = [2.]*89 + [2.]*3 + [1.]*2 + [0.]*(n-94)
assert len(occupation) == n and sum(occupation) == 186
lines += ["#OCC", "* OCCUPATION NUMBERS"]
for first in range(0, n, 5):
    lines.append("".join(f" {x:21.14E}" for x in occupation[first:first+5]))
decoded = []
for line in lines[lines.index("#ORB")+1:lines.index("#OCC")]:
    if line.startswith("*"):
        continue
    assert len(line) % 22 == 0
    decoded.extend(float(line[first+1:first+22]) for first in range(0, len(line), 22))
serialization_error = float(np.max(np.abs(np.asarray(decoded).reshape(n, n).T - converted)))
assert serialization_error < 1e-12, serialization_error
Path("initial.RasOrb").write_text("\n".join(lines) + "\n")
report = dict(normalized_overlap_error=overlap_error, orbital_gram_error=orth_error,
              serialization_error=serialization_error,
              n_basis=n, inactive=89, active=11, active_electrons=8, multiplicity=3,
              ao_permutation=mapping, center_mapping=center_map,
              note="Same occupied frame; Molcas performs its own initial spin-adapted CI solve")
Path("conversion.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report), flush=True)
