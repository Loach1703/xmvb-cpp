"""Same-GUS FeCl2 reference; no new HF guess and no separate CASCI run."""
import json
import os
import time
from pathlib import Path
import numpy as np
import pyscf
from pyscf import gto, lib, mcscf, scf, fci
from guess_reference import read_xmvb_orbitals, complete_orthonormal_mos

assert os.environ.get('SLURM_JOB_ID')
lib.num_threads(int(os.environ['SLURM_CPUS_PER_TASK']))
mol = gto.M(atom='Fe 0 0 0; Cl 0 0 2.136319; Cl 0 0 -2.136319',
            basis='def2-tzvp', cart=True, charge=0, spin=4,
            symmetry=False, verbose=4, max_memory=24000)
assert mol.nao_nr() == 136 and mol.nelectron == 60 and not mol.has_ecp()
s = mol.intor_symmetric('int1e_ovlp')
guess = read_xmvb_orbitals(Path('input.xmi'), 136, 33, np.sqrt(s.diagonal()))
mo = complete_orthonormal_mos(s, guess, 25)
assert np.max(np.abs(mo.T @ s @ mo - np.eye(136))) < 1e-8
np.save('initial_mo.npy', mo)
mc = mcscf.CASSCF(scf.ROHF(mol), 8, (7, 3))
mc.fcisolver = fci.addons.fix_spin_(fci.direct_spin1.FCISolver(mol), shift=1., ss=6.)
mc.fcisolver.conv_tol = 1e-11
mc.conv_tol, mc.conv_tol_grad = 1e-7, 1e-3
mc.max_cycle_macro = 2000
mc.chkfile = 'pyscf.chk'
initial, history = {}, []
native_casci = mc.casci
def record_initial(frame, *args, **kwargs):
    result = native_casci(frame, *args, **kwargs)
    if not initial:
        initial.update(energy=float(result[0]),
            spin_squared=float(fci.spin_op.spin_square(result[2], 8, (7, 3))[0]))
        print('INITIAL', json.dumps(initial), flush=True)
    return result
mc.casci = record_initial
def callback(env):
    row = {'macro': int(env.get('imacro', -1)), 'energy': float(env.get('e_tot', np.nan)),
           'gradient': float(env.get('norm_gorb', np.nan))}
    history.append(row)
    print('MACRO', json.dumps(row), flush=True)
mc.callback = callback
start = time.perf_counter()
mc.kernel(mo)
result = {'converged': bool(mc.converged), 'initial': initial,
          'energy': float(mc.e_tot), 'optimization_wall_s': time.perf_counter()-start,
          'history': history, 'pyscf_version': pyscf.__version__,
          'spin_squared': float(fci.spin_op.spin_square(mc.ci, 8, (7, 3))[0])}
Path('result.json').write_text(json.dumps(result, indent=2))
np.save('final_mo.npy', mc.mo_coeff)
print('FINAL', json.dumps(result), flush=True)
