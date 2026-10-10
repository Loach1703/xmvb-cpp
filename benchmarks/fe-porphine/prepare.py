"""Export an occupied-space-preserving CAS(8e,11o) guess on Slurm.

Fe 3d/4d projectors select two closed and five virtual orbitals. All four
ROHF singly occupied orbitals enter the active space. All runs read the same
frame, including across multiplicities. An optional saved ROHF checkpoint
avoids repeating SCF; it must match this molecule and AO basis.
The primary AO basis is Cartesian, matching xmvb-cpp, not the published
spherical-basis benchmark. No orbitals or spatial irreducible representations
are frozen. See README.md for the comparison's scope.
"""

import hashlib
import json
import os
from pathlib import Path
import sys

assert os.environ.get("SLURM_JOB_ID"), "Run on a Slurm compute node"
import numpy as np
from scipy import linalg
from pyscf import gto, lib, scf
from pyscf.gto.basis import parse_gaussian
from pyscf.tools import molden

lib.num_threads(int(os.environ["SLURM_CPUS_PER_TASK"]))
source = Path(sys.argv[1]).resolve()
atoms = [line.split() for line in Path("geometry.xyz").read_text().splitlines()[2:]]
assert len(atoms) == 37
basis_path = source / "basis/cc-pvdz.gbs"
basis = {a: parse_gaussian.load(str(basis_path), a, optimize=False)
         for a in {row[0] for row in atoms}}
mol = gto.M(atom=[(a[0], [float(x) / 0.529177249 for x in a[1:]]) for a in atoms],
            unit="Bohr", basis=basis, cart=True, spin=4, charge=0,
            symmetry=False, max_memory=120000, verbose=4)
assert mol.nelectron == 186
mf = scf.ROHF(mol)
mf.chkfile = "quintet_rohf.chk"
mf.conv_tol = 1e-10
mf.max_cycle = 200
if len(sys.argv) > 2:
    saved_mol, saved = scf.chkfile.load_scf(sys.argv[2])
    assert saved_mol.cart == mol.cart and saved_mol.spin == mol.spin
    assert saved_mol.charge == mol.charge and saved_mol._basis == mol._basis
    assert saved_mol.atom_charges().tolist() == mol.atom_charges().tolist()
    assert np.max(np.abs(saved_mol.atom_coords() - mol.atom_coords())) < 1e-10
    assert saved["mo_coeff"].shape == (mol.nao_nr(), mol.nao_nr())
    mf.mo_coeff, mf.mo_occ = saved["mo_coeff"], saved["mo_occ"]
    mf.mo_energy, mf.e_tot = saved["mo_energy"], saved["e_tot"]
    mf.converged = True
    import shutil
    shutil.copyfile(sys.argv[2], mf.chkfile)
else:
    mf.kernel()
    if not mf.converged:
        mf = mf.newton()
        mf.max_cycle = 100
        mf.kernel()
assert mf.converged, "Do not export an unconverged ROHF guess"
s = mf.get_ovlp()
spherical = mol.copy()
spherical.cart = False
spherical.build(False, False)
indices = spherical.search_ao_label(["Fe 3d", "Fe 4d"])
assert len(indices) == 10, "Expected two five-dimensional spherical d shells"
# Retain the full Cartesian AO space while using exactly ten d projectors.
imp = mol.cart2sph_coeff()[:, indices]
values, vectors = linalg.eigh(imp.T @ s @ imp)
assert values.min() > 0
imp = imp @ (vectors / np.sqrt(values))
closed = mf.mo_coeff[:, mf.mo_occ == 2]
singly = mf.mo_coeff[:, mf.mo_occ == 1]
virtual = mf.mo_coeff[:, mf.mo_occ == 0]
assert closed.shape[1] == 91 and singly.shape[1] == 4


def select(block, count):
    """Select the largest Fe-projector weights within one occupation sector."""
    projection = imp.T @ s @ block
    weights, rotation = linalg.eigh(projection.T @ projection)
    ordered = block @ rotation[:, ::-1]
    return ordered[:, :count], ordered[:, count:], weights[::-1]


active_closed, core, closed_weights = select(closed, 2)
active_virtual, external, virtual_weights = select(virtual, 5)
n_core = (mol.nelectron - 8) // 2
blocks = [core, np.column_stack((active_closed, singly, active_virtual)), external]
sc = s @ mf.mo_coeff
fock = (sc * mf.mo_energy) @ sc.T
canonical = []
for block in blocks:
    _, rotation = linalg.eigh(block.T @ fock @ block)
    canonical.append(block @ rotation)
mo = np.column_stack(canonical)
error = float(np.max(np.abs(mo.T @ s @ mo - np.eye(mo.shape[1]))))
assert error < 1e-7
occupied = mf.mo_coeff[:, mf.mo_occ > 0]
occupied_leakage = float(occupied.shape[1] - np.sum(
    (mo[:, :n_core + 11].T @ s @ occupied) ** 2))
core_population = float(np.sum(
    (mo[:, :n_core].T @ s @ mf.mo_coeff) ** 2 * mf.mo_occ[None, :]))
assert abs(occupied_leakage) < 1e-7
assert abs(core_population - 2 * n_core) < 1e-7
np.save("initial_mo.npy", mo)
molden.from_mo(mol, "initial.molden", mo)
n_bf = mol.nao_nr()
n_orb = n_core + 11
coeff = mo[:, :n_orb] * np.sqrt(s.diagonal())[:, None]
gus = ["$GUS", " ".join([str(n_bf)] * n_orb)]
for i in range(n_orb):
    gus.append(f"# ORBITAL {i + 1} NAO = {n_bf}")
    for start in range(0, n_bf, 4):
        gus.append(" ".join(f"{coeff[j, i]:.16e} {j + 1}"
                            for j in range(start, min(start + 4, n_bf))))
gus.append("$END")
for name, multiplicity in (("singlet", 1), ("triplet", 3), ("quintet", 5)):
    for iscf in (5, 8, 9):
        lines = [f"Fe porphine {name}; shared occupied-preserving ROHF frame", "$CTRL",
                 f"vbscf basis=cc-pvdz int=libcint ncharge=0 nmul={multiplicity}",
                 "nao=11 nae=8 str=full orbtyp=oeo frgtyp=atom guess=read",
                 f"iscf={iscf} eigensolver=davidson itmax=1000 molden iprint=3",
                 "$END", "$GEO"]
        lines.extend(" ".join(a) for a in atoms)
        lines.extend(["$END", *gus])
        Path(f"{name}-iscf{iscf}.xmi").write_text("\n".join(lines) + "\n")
audit = dict(electrons=186, n_bf=n_bf, inactive_orbitals=n_core,
             active_electrons=8, active_orbitals=11, cartesian=True,
             rohf_energy=float(mf.e_tot), orbital_gram_error=error,
             selected_closed_projector_weights=closed_weights[:2].tolist(),
             selected_virtual_projector_weights=virtual_weights[:5].tolist(),
             core_population=core_population,
             occupied_leakage=occupied_leakage,
             pyscf_version=__import__("pyscf").__version__,
             initial_mo_sha256=hashlib.sha256(Path("initial_mo.npy").read_bytes()).hexdigest(),
             basis_sha256=hashlib.sha256(basis_path.read_bytes()).hexdigest())
Path("initial_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
print("INITIAL_AUDIT", json.dumps(audit), flush=True)
