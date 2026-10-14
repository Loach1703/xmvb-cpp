"""Matched-basis Molcas preparation and AO conversion; Slurm nodes only."""
import json
import os
from pathlib import Path
import re
import sys

assert os.environ.get("SLURM_JOB_ID"), "Run on a Slurm compute node"
import h5py
import numpy as np
from pyscf import scf

root, stage = Path(sys.argv[1]), sys.argv[2]
mol, _ = scf.chkfile.load_scf(str(root / "porphine-reference.chk"))
mol.spin = 0
geometry = []
for line in re.split(r"\$GEO\b", (root / "inputs/singlet/input.xmi").read_text(), flags=re.I)[1].splitlines():
    if line.strip().upper() == "$END":
        break
    fields = line.split()
    if fields:
        geometry.append((fields[0], tuple(float(x)/0.529177249 for x in fields[1:])))
mol.set_geom_(geometry, unit="Bohr", inplace=True)
mol.build(False, False)
assert mol.cart and mol.nao_nr() == 469 and mol.nelectron == 186 and not mol._ecp

if stage == "prepare":
    # Do not substitute the native cc-pVDZ contraction columns. Write the
    # actual reference contractions, combining shells only within the same l.
    basis = ["#Hamiltonian NRH", "#Nucleus UNK", "#Contraction CCC", "#AllElectron AE_"]
    labels = {}
    for symbol in sorted(set(mol.atom_pure_symbol(i) for i in range(mol.natm))):
        atom = next(i for i in range(mol.natm) if mol.atom_pure_symbol(i) == symbol)
        shells = [s for s in range(mol.nbas) if mol.bas_atom(s) == atom]
        blocks = []
        for angular in range(max(mol.bas_angular(s) for s in shells)+1):
            selected = [s for s in shells if mol.bas_angular(s) == angular]
            exponents = sorted(set(float(x) for s in selected for x in mol.bas_exp(s)), reverse=True)
            indices = {exponent: row for row, exponent in enumerate(exponents)}
            columns = []
            for shell in selected:
                for coefficients in mol.bas_ctr_coeff(shell).T:
                    column = np.zeros(len(exponents))
                    for exponent, coefficient in zip(mol.bas_exp(shell), coefficients):
                        column[indices[float(exponent)]] = coefficient
                    columns.append(column)
            blocks.append((angular, exponents, np.column_stack(columns)))
        primitives = "".join(f"{len(exps)}{'spdfgh'[l]}" for l, exps, _ in blocks)
        contractions = "".join(f"{matrix.shape[1]}{'spdfgh'[l]}" for l, _, matrix in blocks)
        label = f"{symbol}.XMVBREF.Reference.{primitives}.{contractions}."
        labels[symbol] = label
        basis += ["/"+label, "Explicit shared reference basis", symbol,
                  f"{mol.atom_charge(atom)}.0 {len(blocks)-1}"]
        for angular, exponents, matrix in blocks:
            basis += [f"* {'spdfgh'[angular]}-type functions", f"{len(exponents)} {matrix.shape[1]}"]
            basis += [f"{x:.15E}" for x in exponents]
            # The native reader uses bounded text records. Six 21-character
            # fields fit in 126 columns without truncating the last contraction.
            basis += ["".join(f"{x:20.12E} " for x in row) for row in matrix]
    Path("XMVBREF").write_text("\n".join(basis)+"\n")
    lines = ["&GATEWAY", "Title = Fe porphine singlet, shared original guess",
             "BASDIR = $CurrDir", "NOCD"]
    for atom in range(mol.natm):
        symbol = mol.atom_pure_symbol(atom)
        xyz = " ".join(f"{x:.16E}" for x in mol.atom_coord(atom))
        lines += ["Basis set", labels[symbol], f"{symbol}{atom+1} {xyz} /bohr",
                  "Cartesian all", "End of basis"]
    gateway = "\n".join(lines)+"\n"
    Path("setup.input").write_text(gateway+"&SEWARD\nONEONLY\n&GUESSORB\n")
    # Keep the same Project/scratch; SEWARD replaces one-electron-only data
    # with conventional four-center integrals before RASSCF consumes them.
    rasscf = """&RASSCF
Title = Fe porphine CAS(8,11) singlet from original shared guess
FileOrb = $CurrDir/initial.RasOrb
Spin = 1
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
    Path("optimize.input").write_text("&SEWARD\n"+rasscf)
    print(json.dumps(dict(n_basis=469, inactive=89, ncas=11, nelecas=8,
        multiplicity=1, cartesian=True, ecp=False, density_fitting=False,
        basis_labels=labels)), flush=True)
    sys.exit(0)
if stage != "convert":
    raise ValueError(stage)

with h5py.File("setup.guessorb.h5", "r") as file:
    assert int(file.attrs["NSYM"]) == 1
    ids = file["BASIS_FUNCTION_IDS"][:]
    n = len(ids)
    sm = file["AO_OVERLAP_MATRIX"][:].reshape(n,n)
    primitive_ids, primitives = file["PRIMITIVE_IDS"][:], file["PRIMITIVES"][:]
    centers, charges = file["CENTER_COORDINATES"][:], file["CENTER_CHARGES"][:]
assert n == 469 and ids.shape == (n,4)
center_map = {}
for center, (xyz, charge) in enumerate(zip(centers, charges),1):
    matches = [i for i in range(mol.natm) if mol.atom_charge(i) == charge
               and np.max(np.abs(mol.atom_coord(i)-xyz)) < 1e-9]
    assert len(matches) == 1, (center,matches)
    center_map[center] = matches[0]
assert sorted(center_map.values()) == list(range(mol.natm))
locations = mol.ao_loc_nr()
components, phases = {}, {}
for center in sorted(set(ids[:,0])):
    for angular in sorted(set(abs(int(row[2])) for row in ids if row[0] == center)):
        for radial in sorted(set(int(row[1]) for row in ids if row[0] == center and abs(row[2]) == angular)):
            rows = primitives[np.all(primitive_ids == [center,angular,radial],axis=1)]
            rows = rows[rows[:,1] != 0]
            rows = rows[np.argsort(-rows[:,0])]
            direction = rows[:,1]/np.linalg.norm(rows[:,1])
            matches = []
            for shell in range(mol.nbas):
                if mol.bas_atom(shell) != center_map[center] or mol.bas_angular(shell) != angular:
                    continue
                for column, coefficients in enumerate(mol.bas_ctr_coeff(shell).T):
                    mask = coefficients != 0
                    exponents = mol.bas_exp(shell)[mask]
                    if len(exponents) != len(rows) or not np.allclose(exponents,rows[:,0],rtol=1e-11,atol=1e-12):
                        continue
                    candidate = coefficients[mask]/np.linalg.norm(coefficients[mask])
                    phase = 1 if candidate@direction > 0 else -1
                    if np.max(np.abs(candidate-phase*direction)) < 1e-9:
                        matches.append((int(locations[shell])+column*(angular+1)*(angular+2)//2,phase))
            assert len(matches) == 1, (center,angular,radial,matches)
            offset,phase = matches[0]
            phases[(center,radial,angular)] = phase
            component = 0
            for x in range(angular,-1,-1):
                for y in range(angular-x,-1,-1):
                    components[(center,radial,angular,x,y,angular-x-y)] = offset+component
                    component += 1
mapping,signs = [],[]
for center,radial,angular,m in ids:
    if angular == 0:
        angular,x,y,z = 0,0,0,0
    elif angular == 1:
        x,y,z = {1:(1,0,0),-1:(0,1,0),0:(0,0,1)}[int(m)]
    else:
        assert angular < 0, "Spherical shell detected"
        angular = -int(angular)
        yz = int(np.floor((np.sqrt(8*(int(m)+angular)+1)-1)/2))
        z = int(m)+angular-yz*(yz+1)//2
        x,y = angular-yz,yz-z
    mapping.append(components[(center,radial,angular,x,y,z)])
    signs.append(phases[(center,radial,angular)])
assert sorted(mapping) == list(range(n))
sp = mol.intor_symmetric("int1e_ovlp")
dp,dm = np.sqrt(sp.diagonal()),np.sqrt(sm.diagonal())
signs = np.asarray(signs)
overlap_error = float(np.max(np.abs(sm/np.outer(dm,dm)-
    (sp/np.outer(dp,dp))[np.ix_(mapping,mapping)]*np.outer(signs,signs))))
assert overlap_error < 1e-8, overlap_error
mo = np.load("shared-mo.npy")
converted = mo[mapping,:]*(signs*dp[mapping]/dm)[:,None]
gram_error = float(np.max(np.abs(converted.T@sm@converted-np.eye(n))))
assert gram_error < 1e-8, gram_error
lines = ["#INPORB 2.2", "#INFO", "* Same original singlet orbitals; matched contractions",
         "       0       1       0",f"{n:8d}",f"{n:8d}","#ORB"]
for orbital in range(n):
    lines.append(f"* ORBITAL    1 {orbital+1:5d}")
    for first in range(0,n,5):
        lines.append("".join(f" {x:21.14E}" for x in converted[first:first+5,orbital]))
occupations = [2.]*89+[2.]*4+[0.]*(n-93)
assert sum(occupations) == 186 and len(occupations) == n
lines += ["#OCC", "* OCCUPATION NUMBERS"]
for first in range(0,n,5):
    lines.append("".join(f" {x:21.14E}" for x in occupations[first:first+5]))
Path("initial.RasOrb").write_text("\n".join(lines)+"\n")
report = dict(normalized_overlap_error=overlap_error, orbital_gram_error=gram_error,
    n_basis=n,inactive=89,ncas=11,nelecas=8,multiplicity=1,ao_permutation=mapping,
    source="Original singlet GUS; no orbital optimization in conversion")
Path("conversion.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report),flush=True)
