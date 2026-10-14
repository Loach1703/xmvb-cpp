# Iron porphine triplet: corrected exact-integral NEO run

Submitted on Hanhai22 on 2026-10-14 as job **1253779**, with 32 CPUs and
192 GiB requested memory. No optimizer or convergence changes were made.

The job uses the same frozen binary as corrected singlet job **1253652**:
`shell-screening-fix-1253607/xmvb-cpp.exe`, SHA-256
`503acdae6d262c4e786c77982401609a33101a30fc531c30c50e317ba31618f5`.
This binary includes the normalized Schwarz shell-screening correction.

- Input is copied unchanged from `inputs/triplet/input.xmi`, including its
  original `$GUS` orbital coefficients and `guess=read`.
- Input SHA-256:
  `79edd76ddc6902f84492592cfd70054a449b403ea150b4429ce5306e415ce8eb`.
- Neutral triplet (`nmul=3`), Cartesian cc-pVDZ, CAS(8,11), 89 inactive
  orbitals, full OEO, full structure space, exact libcint integrals.
- NEO (`iscf=9`) with Davidson and analytic Hessian-diagonal preconditioning.
- Gradient tolerance `1e-3`, energy tolerance `1e-7` Eh, exact-exchange
  memory budget 256 MiB, unchanged from the corrected singlet run.
- Original input iteration limit is 1000; there is no short 4/30-step cap.
- The original-guess PySCF triplet reference energy is
  `-2244.7908349995173` Eh. This reference is already completed, not rerun.

Submission command (from the remote run root):

```sh
sbatch --parsable --job-name=porphine-screened-triplet \
  porphine-screened-converge.sbatch triplet
```

Run directory:
`/home/guqqgroup/xiatao/xmvb-runs/shared-response-20261012/porphine-screened-converge-1253779`.
`output.out` contains the energy/gradient trajectory, `neo.tsv` contains
HVP and keyframe counts, and `resource.time` records wall time and peak RSS
when the computation terminates.
