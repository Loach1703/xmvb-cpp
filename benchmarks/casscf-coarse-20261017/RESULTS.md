# CASSCF coupled coarse correction and recycled initial directions

Parent implementation: `2797d87b`. Exact, single-state CAS(10e,8o), quintet
FeCl2, def2-TZVP Cartesian, 136 AO functions and 25 inactive orbitals.
All measurements used 32 allocated physical cores on Hanhai25 AMD EPYC 9755
nodes. Different nodes were used; these are single native measurements, not
noise-qualified averages. No historical C++ binary was rerun.

## Implementation

- Reuse `mcscf::CoarsePreconditioner`, not a CASSCF-specific inverse.
- One coupled action returns its orbital and CI source images separately.
  The two coarse samples require no additional physical HVP or sigma sweep.
- Keep coarse images at one macro point only. Recycle the last accepted
  search direction across macros; align the relaxed CI eigenvector's sign
  with its finite predictor and project the new initial CI tangent.
- Keep the existing gradient/HVP, orbital retraction, microstep policy,
  rank thresholds and outer energy/gradient tolerances unchanged.
- Cache the small certified Hessian/metric projections; shifted coarse
  applications no longer repeat their full-coordinate projection.

The first integration exposed a numerical-admission defect in the common
coarse model: image QR can amplify source roundoff even when the sample
Gram is still rank-resolved. The projection bound also omitted the existing
full-coordinate dot-product certificate. Candidate samples now satisfy the
existing symmetry certificate for both images before they enter the
optional preconditioner. Uncertifiable samples are not retained. Physical
CIAH actions and their independent symmetry checks remain unchanged. There
is no larger tolerance, additional physical response or fallback solver.

## Native `.xmi` results

These use byte-identical native input decks, the same exact-integral AO
loader and unchanged thresholds (1e-7 Eh and 1e-3 gradient norm).
"Total" is the program's end-to-end wall time; solver time excludes input
and AO integral generation. RSS is GNU time's process high-water mark.

| Version | Steps | HVP | Final energy / Eh | Solver / s | Total / s | RSS / MiB |
|---|---:|---:|---:|---:|---:|---:|
| Before, 291886 | 25 | 355 | -2181.624464652800 | 4.80796 | 5.10977 | 400.46 |
| After, 291946 | 18 | 255 | -2181.624464862390 | 3.50281 | 3.79407 | 407.47 |

Steps decrease 28%; HVP decreases 28.2%; end-to-end wall decreases 25.7%.
The retained coarse model adds about 7 MiB to the measured process peak,
not a per-thread copy of AO integrals or large matrices.

## Identical AO/MO fixture and correctly threaded PySCF

Both solvers read the same numerical integral arrays and full initial frame,
SHA-256 `801c0fddf648e281e0453fe00a0813a525514b98ff047908fbaad36594232190`.
Initial CASCI energy is -2181.617631912353 Eh (C++ difference 4.55e-13 Eh).
The fixture has unscreened AO integrals and a larger C++ AO graph than the
native loader. Its times/RSS must be compared within this path, not mixed
with the native `.xmi` row above.

| Method | Steps | HVP/action | Final energy / Eh | Solver / s | Process / s | RSS / MiB |
|---|---:|---:|---:|---:|---:|---:|
| C++ before, 291883 | 20 | 277 | -2181.624464828963 | 12.60521 | 13.30 | 1664.00 |
| C++ after, 291944 | 16 | 209 | -2181.624464858833 | 12.10532 | 12.81 | 1664.00 |
| PySCF standard AH, 291935 | 34 | not exposed | -2181.624464771845 | 14.95014 | 15.70 | 776.25 |
| PySCF coupled Newton, 291936 | 24 | 324 | -2181.624464841431 | 24.74374 | 27.01 | 776.54 |

The exact-fixture C++ HVP count falls 24.5%, but solver wall only falls 4.0%.
This round does not establish where the remaining time is spent; fewer HVPs
must not be presented as an equal fractional wall-time improvement. The
fixture C++ peak is unchanged but remains larger than PySCF's packed-ERI
representation. The final C++ energy gaps are 8.70e-8 Eh to standard AH and
1.74e-8 Eh to coupled Newton. All reported final roots have the intended
quintet spin; PySCF reports S(S+1)=6.

The original PySCF 291883/291887 timings were invalid parallel baselines:
their Python worker threads all inherited one CPU. The corrected jobs keep
Slurm's 32-core affinity but disable Python-internal OpenMP pinning. Matched
one-macro evidence and the exact environment are in
`benchmarks/casscf-threads-20261017/RESULTS.md`. The 34 AH macros here are not
its 117 callback records. AH does not expose a cumulative HVP count through
the current callback, so no number is invented.

## Qualification

Build 291943 passed all 12 selected tests, including coupled split-action
equivalence, pure-spin tangent projection, analytic directional derivatives,
finite keyframe pullbacks, and the common balancing inverse. New tests check
full-coordinate roundoff admission, image/metric rejection and recovery,
without modifying the source physical action.

Final CASSCF F2 smoke 291953 converged in 3 steps / 14 HVP to
-198.761111541649 Eh, solver 0.14367 s and 56 MiB RSS; the stored pre-change
native run used 4 steps / 18 HVP. Shared-module VB regressions 291948/291950
still converge FeCl2 in 13 steps (-2181.624464871767 Eh) and F2 in 3 steps
(-198.7511558305 Eh). Their HVP totals are not exposed in these plain logs.

Two F2 array attempts were cancelled by the scheduler before application
startup on anode018; the final retry on anode029 completed. Failed initial
coarse integrations and scalar diagnostic runs are retained, not counted as
performance results. `results.csv` is generated by the common `collect.py`;
artifact IDs distinguish development attempts and the final qualified run.

Local raw root: `results/diagnostics/casscf-coarse-20261017/` (ignored).
Remote root:
`/home/guqqgroup/taoxia/work/xmvb-ciah-subspace-20261017/benchmarks/casscf-coarse-20261017/`.
