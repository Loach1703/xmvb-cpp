"""Prepare a small hydrogen-chain orbital map; no SCF optimization."""
import itertools
import os
from pathlib import Path
import sys
assert os.environ.get('SLURM_JOB_ID')
import numpy as np
from pyscf import gto, lib
lib.num_threads(1)
n = int(sys.argv[1])
mol = gto.M(atom=[('H',(0,0,1.4*i)) for i in range(n)],
            basis='sto-6g', spin=0, verbose=0)
r = np.linalg.cholesky(mol.intor('int1e_ovlp')).T
with Path('orbitals.txt').open('w') as handle:
    handle.write(str(n)+'\n')
    np.savetxt(handle, r, fmt='%.17g')
if n == 12:
    strings = list(itertools.combinations(range(n),n//2))
    rows = np.asarray(strings)
    alpha = [tuple(2*g+bit for g,bit in enumerate(bits))
             for bits in itertools.product((0,1),repeat=n//2)]
    beta = [tuple(sorted(set(range(n))-set(a))) for a in alpha]
    def exterior(columns):
        return np.column_stack([np.linalg.det(
            r[rows[:,:,None],np.asarray(c)[None,None,:]]) for c in columns])
    y = (exterior(alpha)*2.**(-n/4))@exterior(beta).T
    masks = np.asarray([sum(1<<p for p in a) for a in strings],dtype=np.uint64)
    keys = (masks[:,None] | (masks[None,:]<<np.uint64(n))).ravel()
    values = y.ravel()
    order = np.argsort(keys)
    reference = np.empty(len(order),dtype=[('key','<u8'),('value','<f8')])
    reference['key'],reference['value'] = keys[order],values[order]
    reference.tofile('reference.bin')
