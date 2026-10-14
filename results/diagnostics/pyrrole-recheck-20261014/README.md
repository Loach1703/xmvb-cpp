# Pyrrole dimer 3.5 Angstrom recheck

User attachments are preserved in `inputs/oeo.xmi` and `inputs/blo.xmi`.
Both specify CAS(11e,10o), doublet cation, two equally weighted states,
aug-cc-pVDZ, exact Libcint integrals and read-in orbitals. BLO uses HAO
constraints; OEO uses full AO coefficients. Do not require constrained BLO
to match unconstrained CASSCF energies.

Hanhai22 array **1255187** runs four independent 32-core tasks:

| Index | Orbitals | Optimizer |
|---|---|---|
| 0 | OEO | input ISCF=8 block-LBFGS |
| 1 | BLO | input ISCF=8 block-LBFGS |
| 2 | OEO | CLI override NEO |
| 3 | BLO | CLI override NEO |

Binary: shared-response-20261012/horizontal-build-1255164/xmvb-cpp.exe,
SHA256 `bc2003c898dddf977d76ac6bb6e4a5fd03f0700dd9957595c2eb73e830cfcfbb`.
This is the CIAH interior-on-relaxed-rejection implementation archived in
commit 154ba75c. No new code change or threshold tuning in this recheck.

Existing PySCF historical same-guess reference is recorded in
`../pyrrole-molcas-20261002/pyscf-same-guess/result.json`:
initial average -417.49257149477376 Eh; final roots
-417.51528410817974 and -417.48042370161835 Eh; final average
-417.49785390489905 Eh. Its recorded tolerances are stricter and its timing
is from Hanhai25, so do not present it as a matched timing benchmark.
The new OEO attachment matches the historical OEO block-LBFGS CTRL and GEO
sections, but its GUS section is not textually identical. Therefore the old
PySCF run is an energy reference, not yet a verified same-guess comparison
for this new attachment.

Remote artifacts: `/home/guqqgroup/xiatao/xmvb-runs/pyrrole-recheck-20261014/run-1255187/`.
All four tasks were RUNNING at the first check. No convergence result yet.
