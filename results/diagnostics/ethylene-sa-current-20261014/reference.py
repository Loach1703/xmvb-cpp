"""Same-orbital-span Cartesian full-CAS reference, ground or equal-weight SA."""
import json
import os
from pathlib import Path
import re
import sys
import time

import numpy as np
from pyscf import gto, lib, mcscf, scf
from pyscf.tools import molden
from FeCl2_pyscf import complete_orthonormal_mos, read_xmvb_orbitals

if not os.environ.get('SLURM_JOB_ID'):
    raise SystemExit('Run on a Slurm compute node only')
lib.num_threads(int(os.environ['OMP_NUM_THREADS']))
def record_affinity(stage):
    """Record real native-thread masks, not just the Slurm CPU request."""
    masks = {}
    for task in Path('/proc/self/task').iterdir():
        try:
            status = (task / 'status').read_text()
        except FileNotFoundError:
            continue
        mask = re.search(r'^Cpus_allowed_list:\s*(.+)$', status, re.M)[1]
        masks[mask] = masks.get(mask, 0) + 1
    payload = {'stage': stage, 'threads_by_cpu_mask': masks,
               'pyscf_threads': lib.num_threads()}
    Path('affinity-' + stage + '.json').write_text(json.dumps(payload, indent=2))
    print('AFFINITY ' + json.dumps(payload), flush=True)
    if sum(masks.values()) > 1 and len(masks) == 1:
        mask = next(iter(masks))
        if ',' not in mask and '-' not in mask and lib.num_threads() > 1:
            raise RuntimeError('All native threads bound to one CPU')
record_affinity('imports')
source = Path(sys.argv[1])
text = source.read_text()
def integer(key, default=None):
    match = re.search(r'\b' + key + r'\s*=\s*(\d+)', text, re.I)
    if match:
        return int(match[1])
    if default is None:
        raise ValueError(f'Missing {key}')
    return default

assert re.search(r'\bstr\s*=\s*full\b', text, re.I)
assert re.search(r'\borbtyp\s*=\s*oeo\b', text, re.I)
geometry = re.search(r'\$geo\s*\n(.*?)\$end', text, re.I | re.S)[1]
atoms = [line.split('#')[0].strip() for line in geometry.splitlines() if line.split('#')[0].strip()]
basis = re.search(r'\bbasis\s*=\s*([\w+-]+)', text, re.I)[1]
spin = integer('nmul', 1) - 1
n_active = integer('nao')
n_electron = integer('nae')
weights_match = re.search(r'wstate\(1\)\s*=\s*([^\n]+)', text, re.I)
weights = [float(x) for x in weights_match[1].strip().split(',')] if weights_match else [1.]
weights = np.array(weights) / sum(weights)
# Read only the basis/shell layout, never the final orbital columns. Diffuse
# shells occur after the parent basis in C++; named PySCF bases sort by l.
layout = Path('../cpp/input.molden')
if not layout.exists():
    raise ValueError('Missing C++ basis layout: cannot certify the AO ordering')
mol = molden.load(str(layout))[0]
mol.charge = integer('ncharge', 0)
mol.spin = spin
mol.symmetry = False
# Restore the exact input coordinates (Molden writes rounded coordinates).
# The ndarray set_geom_ path updates only coordinates, preserving shell order.
coordinates = np.array([[float(x) for x in atom.split()[1:4]] for atom in atoms])
mol.set_geom_(coordinates / 0.529177249, unit='Bohr', symmetry=False, inplace=True)
mol.verbose = 4
mol.max_memory = 48000
assert mol.cart
n_core = (mol.nelectron - n_electron) // 2
assert 2 * n_core + n_electron == mol.nelectron
overlap = mol.intor_symmetric('int1e_ovlp')
guess = read_xmvb_orbitals(source, mol.nao, n_core + n_active, np.sqrt(np.diag(overlap)))
mo = complete_orthonormal_mos(overlap, guess, n_core)
mf = scf.ROHF(mol)
mf.max_memory = 48000
electrons = ((n_electron + spin)//2, (n_electron - spin)//2)
def configure(model):
    model.fcisolver.spin = spin
    model.fix_spin_(ss=spin/2*(spin/2+1))
    if len(weights) > 1:
        model = model.state_average_(weights.tolist())
    return model

initial = configure(mcscf.CASCI(mf, n_active, electrons))
initial_energy = float(initial.kernel(mo)[0])
record_affinity('initial')
print(f'INITIAL_ENERGY = {initial_energy:.15f}', flush=True)
cpp_output = Path('../cpp/output.out').read_text()
iteration_table = cpp_output.split('DELTA_E')[1]
cpp_initial = float(re.search(r'^\s*0\s+([-+\d.]+)', iteration_table, re.M)[1])
initial_difference = initial_energy - cpp_initial
print(f'INITIAL_ENERGY_DIFFERENCE = {initial_difference:.12e}', flush=True)
if abs(initial_difference) > 1e-6:
    raise ValueError('Initial energies disagree: do not run an invalid comparison')
driver = sys.argv[2] if len(sys.argv) > 2 else 'newton'
model = configure(mcscf.CASSCF(mf, n_active, electrons))
if driver == 'newton':
    model = model.newton()
elif driver != 'ah':
    raise ValueError(f'Unknown PySCF driver: {driver}')
model.conv_tol = 1e-7
model.conv_tol_grad = 1e-3
model.max_cycle_macro = 2000
history = []
def callback(env):
    row = {'macro': int(env.get('imacro', len(history)+1)),
           'energy': float(env.get('e_tot', np.nan)),
           'gradient': float(env.get('norm_gall', env.get('norm_gorb', np.nan)))}
    history.append(row)
    print('MACRO ' + json.dumps(row), flush=True)
model.callback = callback
start = time.perf_counter()
energy = float(model.kernel(mo)[0])
elapsed = time.perf_counter() - start
record_affinity('final')
result = {'initial_energy': initial_energy, 'energy': energy,
          'driver': driver,
          'converged': bool(model.converged), 'optimization_seconds': elapsed,
          'macro': history[-1]['macro'] if history else None,
          'basis': basis, 'n_bf': mol.nao, 'n_active': n_active,
          'n_electron': n_electron, 'n_core': n_core, 'spin': spin,
          'weights': weights.tolist(), 'state_energies': np.asarray(model.e_states if len(weights)>1 else [energy]).tolist(),
          'history': history}
Path('reference.json').write_text(json.dumps(result, indent=2))
np.save('final_mo.npy', model.mo_coeff)
print('RESULT ' + json.dumps(result), flush=True)
