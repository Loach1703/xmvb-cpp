"""Export an approved core/active/virtual MO frame to both solvers.

No HF, AVAS, CASCI or optimization is run here. Execute inside Slurm only.
The checkpoint and MO file must come from one approved preparation run.
All basis functions must use the supplied Gaussian basis file and Cartesian AOs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("orbitals", type=Path)
    parser.add_argument("geometry", type=Path)
    parser.add_argument("basis", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--ncas", type=int, required=True)
    parser.add_argument("--nelecas", type=int, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Orbital validation must run in a Slurm allocation")
    import numpy as np
    from pyscf import gto, lib, scf
    from pyscf.gto.basis import parse_gaussian
    lib.num_threads(int(os.environ["SLURM_CPUS_PER_TASK"]))
    mol, _ = scf.chkfile.load_scf(str(args.checkpoint))
    assert mol.cart and not mol.has_ecp() and not mol.symmetry
    atoms = [line.split() for line in args.geometry.read_text().splitlines()[2:] if line.strip()]
    assert len(atoms) == int(args.geometry.read_text().splitlines()[0]) == mol.natm
    basis = {a[0]: parse_gaussian.load(str(args.basis), a[0], optimize=False) for a in atoms}
    rebuilt = gto.M(atom=[(a[0], [float(x)/0.529177249 for x in a[1:]]) for a in atoms],
                    unit="Bohr", basis=basis, cart=True, charge=mol.charge,
                    spin=mol.spin, symmetry=False, verbose=0)
    assert rebuilt._basis == mol._basis, "Basis differs from checkpoint"
    assert np.max(np.abs(rebuilt.atom_coords()-mol.atom_coords())) < 1e-10
    assert all(rebuilt.atom_charge(i) == mol.atom_charge(i) for i in range(mol.natm))
    n_bf = mol.nao_nr()
    mo = np.load(args.orbitals, allow_pickle=False)
    assert mo.shape == (n_bf, n_bf) and np.isrealobj(mo) and np.isfinite(mo).all()
    assert (mol.nelectron-args.nelecas) % 2 == 0
    n_core = (mol.nelectron-args.nelecas)//2
    assert n_core >= 0 and 0 < args.ncas <= n_bf-n_core
    assert (args.nelecas-mol.spin) % 2 == 0
    assert 0 <= (args.nelecas-mol.spin)//2 <= (args.nelecas+mol.spin)//2 <= args.ncas
    s = mol.intor_symmetric("int1e_ovlp")
    gram_error = float(np.max(np.abs(mo.T@s@mo-np.eye(n_bf))))
    assert gram_error < 1e-8
    n_orb = n_core+args.ncas
    coeff = mo[:, :n_orb]*np.sqrt(s.diagonal())[:, None]
    lines = ["Shared approved initial orbitals", "$CTRL",
             f"vbscf basis=shared int=libcint ncharge={mol.charge} nmul={mol.spin+1}",
             f"nao={args.ncas} nae={args.nelecas} str=full orbtyp=oeo frgtyp=atom guess=read",
             "iscf=9 eigensolver=davidson itmax=1000 molden iprint=3", "$END", "$GEO"]
    lines += [" ".join(a) for a in atoms]+["$END", "$GUS", " ".join([str(n_bf)]*n_orb)]
    for i in range(n_orb):
        lines.append(f"# ORBITAL {i+1} NAO = {n_bf}")
        for start in range(0,n_bf,4):
            lines.append(" ".join(f"{coeff[j,i]:.17f} {j+1}" for j in range(start,min(start+4,n_bf))))
    lines.append("$END")
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/"input.xmi").write_text("\n".join(lines)+"\n")
    shutil.copyfile(args.basis,args.output/"shared.gbs")
    shutil.copyfile(args.checkpoint,args.output/"guess.chk")
    shutil.copyfile(args.orbitals,args.output/"initial_mo.npy")
    data = dict(ncas=args.ncas, nelecas=args.nelecas, ncore=n_core, spin=mol.spin,
                gram_error=gram_error, symmetry=False, cartesian=True,
                frame_order="inactive, active, virtual", selection="externally approved; not inferred here")
    data["sha256"] = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.output.iterdir()}
    (args.output/"shared.json").write_text(json.dumps(data,indent=2)+"\n")


if __name__ == "__main__":
    main()
