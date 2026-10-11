# NEO projected-model reuse and fine-grained CI scheduling

## Scope and correctness

Baseline production source: `5ace1f4` (the later `d869641` commit only adds
benchmark documentation). Both executables use exact Libcint, NEO/Davidson,
the same MnF2 input, and unchanged gradient/energy tolerances of 1e-3 and
1e-7 hartree. No Hessian, KKT, symmetry, range, or metric-inverse tolerance
was weakened. No integral or physical-matrix replica was added per thread.

The accepted-point workspace now owns the small projected factors:

| Cached quantity | Invalidated by |
| --- | --- |
| Gradient physical dual norm | A new accepted-point workspace |
| Orbital metric whitening | Appending an orbital basis column |
| Structure pseudoinverse and null vectors | Appending a structure basis column |
| Relaxed projected spectrum and spectral gradient | Appending either basis |

Basis columns are append-only within the workspace, so their dimensions are
valid revision keys. Retained old projected entries never change when a new
column is appended. A radius-only retry reuses all factors, but still solves
its new scalar trust problem and certifies its new full coupled KKT residual.
Symmetry and adjoint checks remain mandatory when their inputs change.

Optional lowest-Ritz residual contractions are deferred during a regular
KKT refinement. They are evaluated on every final candidate and on every
iteration requesting a curvature certificate, including hard cases. The
same final diagnostics and stopping requirements are retained. The gradient
dual norm is reused because both the gradient and physical metric are fixed
within this workspace; candidate residual norms are not substituted or skipped.

The remaining complete projected spectrum uses `LAPACKE_dsyevd` in the owned
projected matrix. The configured Eigen backend previously called `dsyev`.
This removes the QR full-spectrum path and its subsequent eigenvector copy;
it does not approximate eigenpairs or change spectral/trust-region criteria.

Short CI connection/generator loops use either serial execution or the
existing requested team. Parallel execution requires at least one team-width
of items per worker, using an overflow-safe size comparison. This is a
general granularity policy, not an accuracy condition or a molecule-specific
cutoff. Heavy sigma and tiled integral-adjoint GEMMs keep their existing
parallel planner; their cost cannot be inferred from row count alone.

## Matched Hanhai25 results

Job `281273`, node `anode020`, one 32-core Slurm allocation. Current and
baseline binaries run consecutively for each of 32, 16, and 8 threads, with
identical compact affinity and single-threaded BLAS. Each configuration is
measured once. Wall time is the application's end-to-end timer, including
input/integral preparation. Peak RSS is GNU time's process maximum in MiB.

| Cores | Version | Wall, s | Steps | Orbital HVPs | Energy, hartree | Peak RSS, MiB | Average CPU cores |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | Baseline | 26.430535 | 14 | 669 | -1348.893353119312 | 406.1 | 7.62 |
| 8 | Current | 26.144084 | 14 | 669 | -1348.893353119314 | 407.6 | 7.76 |
| 16 | Baseline | 17.234356 | 14 | 669 | -1348.893353119316 | 460.8 | 14.77 |
| 16 | Current | 16.626610 | 14 | 669 | -1348.893353119309 | 473.8 | 15.21 |
| 32 | Baseline | 13.572954 | 14 | 669 | -1348.893353119316 | 495.5 | 28.49 |
| 32 | Current | 12.682297 | 14 | 669 | -1348.893353119323 | 486.9 | 29.47 |

Time reductions are 1.1%, 3.5%, and 6.6%. Energy differences are below
1e-11 hartree. This is a modest total-time improvement, not an order-of-magnitude
gain. RSS is not uniformly lower: changes are +1.5, +12.9, and -8.7 MiB.
Average CPU includes spinning, not just useful work. The 32-core CPU set spans
two NUMA domains; 8/16 cores span one. The node is not exclusively allocated.

| Stage, seconds | Baseline 8 | Current 8 | Baseline 16 | Current 16 | Baseline 32 | Current 32 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Projection factor/model solve | 1.310 | 0.600 | 1.352 | 0.612 | 1.436 | 0.652 |
| Projection residual/certification | 1.465 | 1.467 | 1.477 | 1.448 | 1.758 | 1.777 |
| Physical metric solve | 1.267 | 1.281 | 1.310 | 1.303 | 1.375 | 1.387 |
| CI generator | 0.00755 | 0.00775 | 0.01337 | 0.02163 | 0.12867 | 0.02194 |
| CI connection construction | 0.00030 | 0.00060 | 0.01149 | 0.00064 | 0.07305 | 0.00070 |

Scopes are nested; do not sum these rows. Factor time drops by about 54-55%.
Residual certification time is essentially unchanged. At 32 cores the
physical metric solve accounts for about 1.39 s within that residual stage;
this has not been solved by merely removing repeated calculations.

There are 698 projected iterations, but only 669 metric factorizations,
39 structure factorizations, and 694 relaxed-spectrum solves in the current
run. Most spectrum solves are required because the orbital subspace genuinely
changes. Reuse alone cannot remove them. Metric inverse calls fall from 756
to 752: four duplicate accepted-gradient norm evaluations are removed, not
hundreds of genuinely distinct residual certificates. The spectral backend
change is the principal measured factor-stage gain.

## Rejected implementation

Job `281269` tested varying intermediate CI team sizes. It remained correct
(14 steps, 669 HVPs), but its 32-core wall time regressed from 13.61 to 30.83 s;
CI-generator time increased to 2.33 s. That policy was removed, not retained
as an alternate path. Job `281270` tested the fixed-team/serial policy and
reuse with the old QR spectrum; its 32-core wall time was 13.26 s. Final job
`281273` adds the divide-and-conquer spectrum. These developmental runs are
not additional statistical repetitions of the final configuration.

## Verification and artifacts

All six regression executables passed at 1, 8, and 32 threads under Slurm:
`test_active_two_electron_sparse`, `test_ri_ao_h1e_hvp`,
`test_orthogonal_direct_ci`, `test_neo_response_solver`, `test_neo`, and
`test_neo_globalization`. Tests include dense Schur/KKT references, singular
structure blocks, trust-radius workspace reuse, and final deferred-curvature
diagnostics. The radius-reuse test additionally counts accepted-gradient
dual-norm evaluations and requires exactly one within the workspace.

Reproduce with `mnf2-storage.sbatch SOURCE_DIR BASELINE_EXE reverse`.
The helper saves build flags, input/binary hashes, allocation, source hashes,
regression logs, phase timings, affinity, and GNU-time resources.
Local source hashes match the compute job's saved source hashes.

Raw evidence: `results/studies/locality-20261011/mnf2-storage-281273/`;
parsed evidence: `neo-reuse-281273.csv` and `neo-reuse-phases-281273.csv`.
Remote root: `/home/guqqgroup/taoxia/xmvb-runs/native-20261011/`.

Baseline executable SHA256:
`1c2fb68686f7333158aba37e21a12a2cf3243f353cb003cc1ce8ff40098d4f2b`.
Current executable SHA256:
`c76c5401d4faf53a71d08e83a33ceca5c77763b196e7cab526dccb15e3cde06f`.
Input SHA256:
`79034c0fdb5baee2432e81d90c7316de2b0bf2edfe35166b184f2c28f761a822`.

Remaining priorities are the strict physical metric inverse and the essential
changing-subspace spectrum. Optimize their actual structure or data movement;
do not claim that their remaining calls are redundant or loosen certification.
