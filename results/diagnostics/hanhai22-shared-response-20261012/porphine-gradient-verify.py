"""Compare full orbital gradients for the same corrected C++ wavefunction."""
import json
import os
from pathlib import Path
import re
import sys
import numpy as np
from pyscf import fci, lib, mcscf, scf

assert os.environ.get("SLURM_JOB_ID")
root, cpp = map(Path,sys.argv[1:3])
lib.num_threads(int(os.environ["SLURM_CPUS_PER_TASK"]))
mol,saved = scf.chkfile.load_scf(str(root/"porphine-reference.chk"))
mol.spin = 0
geometry = []
for line in re.split(r"\$GEO\b",(cpp/"endpoint.xmi").read_text(),flags=re.I)[1].splitlines():
    if line.strip().upper() == "$END":
        break
    fields = line.split()
    if fields:
        geometry.append((fields[0],tuple(float(x)/0.529177249 for x in fields[1:])))
mol.set_geom_(geometry,unit="Bohr",inplace=True)
mol.build(False,False)
assert mol.cart and mol.nao_nr() == 469 and mol.nelectron == 186 and not mol._ecp
mf = scf.ROHF(mol)
mf.mo_coeff,mf.mo_occ = saved["mo_coeff"],saved["mo_occ"]
mf.mo_energy,mf.e_tot,mf.converged = saved["mo_energy"],saved["e_tot"],True
s = mf.get_ovlp()
d = np.sqrt(s.diagonal())
mo = np.load(root/"porphine-terminal-1253506/pyscf/initial_mo.npy")
prefix = cpp/"audit_orbital_gradient_"
raw = np.fromfile(str(prefix)+"accepted_raw_f64.bin").reshape((469,100),order="F")
previous = np.fromfile(root/"porphine-terminal-1253506/cpp/audit_orbital_gradient_accepted_raw_f64.bin").reshape((469,100),order="F")
assert np.max(np.abs(raw-previous)) < 1e-13, "The saved MO frame is not this raw point"
gradient = np.fromfile(str(prefix)+"packed_gradient_f64.bin").reshape((469,100),order="F")
def masks(side):
    return [sum(1<<int(i) for i in line.split()) for line in
        Path(str(prefix)+side+"_strings.txt").read_text().splitlines()]
alpha,beta = masks("alpha"),masks("beta")
ia,ib = {m:i for i,m in enumerate(alpha)},{m:i for i,m in enumerate(beta)}
strings = list(map(int,fci.cistring.make_strings(range(11),4)))
carrier = np.fromfile(str(prefix)+"carrier_f64.bin").reshape((330,330),order="F")
ci = carrier[np.ix_([ia[m] for m in strings],[ib[m] for m in strings])]
ci /= np.linalg.norm(ci)
mc = mcscf.CASSCF(mf,11,(4,4))
mc.max_memory = 120000
mc.fcisolver.spin = 0
mc.fix_spin_(ss=0)
k = mo[:,:100].T@s@(raw/d[:,None])
representation_error = float(np.max(np.abs(mo[:,:100]@k-raw/d[:,None])))
a = mo.T@(d[:,None]*gradient)@k.T
image = np.zeros((469,469))
image[:,:100] = a
cpp_g = mc.pack_uniq_var(image-image.T)
py_g = 2*mc.get_grad(mo,fci.direct_spin1.make_rdm12(ci,11,(4,4)))
np.save("cpp_rotation_gradient.npy",cpp_g)
np.save("pyscf_rotation_gradient.npy",py_g)
error = mc.unpack_uniq_var(cpp_g-py_g)
report = dict(representation_error=representation_error,
    cpp_norm=float(np.linalg.norm(cpp_g)),pyscf_norm=float(np.linalg.norm(py_g)),
    error_norm=float(np.linalg.norm(cpp_g-py_g)),error_inf=float(np.max(np.abs(cpp_g-py_g))),
    relative_error=float(np.linalg.norm(cpp_g-py_g)/np.linalg.norm(py_g)),
    core_active_error=float(np.linalg.norm(error[89:100,:89])),
    core_virtual_error=float(np.linalg.norm(error[100:,:89])),
    active_virtual_error=float(np.linalg.norm(error[100:,89:100])),
    spin_square=float(fci.spin_op.spin_square(ci,11,(4,4))[0]))
Path("comparison.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2),flush=True)
