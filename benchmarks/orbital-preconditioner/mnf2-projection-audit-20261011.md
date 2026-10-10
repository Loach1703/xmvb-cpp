# MnF2: NEO projection attribution and rejected optimizations

## Scope

Single-node OpenMP, exact Libcint integrals, NEO/Davidson; no MPI or RI
changes. Source baseline: 49a2d8f. Native runs execute under Hanhai25 Slurm,
with compact allocated-core affinity and single-thread BLAS. Compilation uses
the same Release flags as the preceding shared-link audit. No accuracy,
convergence, subspace, or molecule-specific thread thresholds are changed.

The preceding gprofng capture (280820) identifies the exact AO-pair and AO-H1E
loops and OpenMP barriers as investigation targets. Its timer-period warning
precludes treating sample percentages as accurate wall-time shares. This
round uses native scoped timing to test specific hypotheses.

## The projected-solve timer was not an eigensolver timer

New opt-in phases separate projected factorization/model solution from full
residual reconstruction and certification. The latter includes the physical
dual-norm metric solve and, when needed, final structure-response correction.
The factor phase includes the scalar trust-region solver. These phases are
nested: neither metric nor secular time should be added to their parent totals.

Job 280825, anode009, unchanged numerical implementation:

| Threads | Projected solve / s | Factor/model / s | Residual certification / s | Secular solve / s |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 2.531 | 1.176 | 1.355 | 0.00683 |
| 16 | 2.636 | 1.246 | 1.389 | 0.00699 |
| 32 | 3.204 | 1.407 | 1.797 | 0.00812 |

There are 698 projected solves and 756 metric solves. The corresponding total
metric-solve times are 1.120, 1.196, and 1.344 seconds. These are not all confined
to projected residual certification. Trust-region scalar root solving is not
a useful optimization target for this workload.

## Structure-inverse reuse: no demonstrated factorization benefit

The structure projected block is unchanged when only the orbital basis expands
or the radius changes. A trial retained its spectral inverse/null vectors until
the append-only structure basis grew; all adjoint/range/KKT checks remained.

Job 280849, anode009; paired native runs with identical instrumentation and
CPU sets, candidate first at each thread count:

| Threads | Baseline wall / s | Trial wall / s | Baseline RSS / MiB | Trial RSS / MiB |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 27.078 | 26.554 | 442.6 | 410.4 |
| 16 | 17.936 | 17.800 | 500.8 | 469.4 |
| 32 | 14.454 | 14.284 | 508.4 | 497.2 |

At 16 threads, the targeted factor phase is 1.2458 versus 1.2497 seconds:
effectively unchanged. Small whole-program/RSS differences do not establish
that this cache improved the target. The cache implementation was removed.
NEO, coupled-response and spectral trust-region unit tests passed at 1/8/32
threads in job 280850, including workspace reuse after radius changes and
singular structure-block range checks.

## SIMD accumulator experiment

A second trial used eight local double accumulators, matching one AVX-512 or
two AVX2 vectors, keeping accumulation local over a complete sparse row. It
revisited CSR metadata for each channel block rather than writing the image
after each four-edge group. Physical arrays stayed shared and the per-thread
scratch was only 64 bytes. Nevertheless, this tradeoff regressed on MnF2:
less output traffic did not compensate for the repeated sparse traversal.

The trial was removed, not retained as a second path or fallback. The original
CSR SIMD contraction and all preceding shared-link/parallel improvements
remain. Precise paired results and resources are saved under job 280851.

| Threads | Baseline wall / s | Trial wall / s | AO-pair baseline / s | AO-pair trial / s | Baseline RSS / MiB | Trial RSS / MiB |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 27.418 | 40.356 | 8.555 | 21.910 | 404.0 | 442.8 |
| 16 | 17.669 | 24.584 | 4.682 | 11.971 | 488.1 | 453.2 |
| 32 | 14.912 | 18.022 | 3.079 | 7.153 | 495.0 | 494.4 |

Baseline mean accounted CPU cores are 7.65, 14.89 and 28.89; trial values are
7.77, 15.22 and 29.49. Increased CPU accounting does not imply more useful
work. Relative to each version's 8-thread result, 16/32-thread scaling
efficiency is 77.6%/46.0% for baseline versus 82.1%/56.0% for trial. The trial
scales better relative to its own slower 8-thread run, yet is slower at every
core count; scaling efficiency alone is not a performance objective. These
are 8-thread-reference efficiencies, not single-thread-reference efficiencies.

## Retained changes and evidence

Only phase subdivision and a compute-node NEO verification script remain in
production/reusable tooling. Neither modifies the optimization mathematics.
Job 280854 rebuilds the restored production executable and passes 18 test
executions: NEO, coupled-response, spectral trust region, exact sparse active
integrals, RI AO-H1E primitives and orthogonal direct CI, each at 1/8/32 threads.
Every completed native run above converges in 14 iterations with 669 orbital
HVPs; job 280849 energies span only about 8e-12 hartree around
-1348.89335311931 hartree. This round does not claim a production speedup.

Raw local logs, binary/input hashes and parsed CSVs:
`results/studies/locality-20261011/mnf2-balance-280825/`,
`mnf2-storage-280849/`, `neo-verify-280850/`, `mnf2-storage-280851/`,
and `neo-verify-280854/`.
Remote root: `/home/guqqgroup/taoxia/xmvb-runs/native-20261011/`.
The immutable trial/baseline executables remain in their remote job directories.

Next investigate the shared exact AO sparse traversal and its consumers using
these phase measurements. Prefer eliminating traversal/data movement over
adding an unproven cache or increasing the number of small parallel regions.
