# Analytic physical-metric blocks and HVP parallel-region audit

## Method and accuracy

Baseline: 7388fc6. Exact/libcint, Davidson, NEO, MnF2, 32 OpenMP threads;
gradient/energy criteria remain 1e-3 / 1e-7 Eh. Builds, unit tests, native
benchmarks, sampling and report generation ran only in Hanhai25 Slurm jobs.

Section 15.1 of [the theory note](../../article/matrix_free_vbscf_orbital_optimization_theory.md)
derives the exact orbital-diagonal blocks of the coupled physical metric for
both individual-ray and complete-active-span charts. These are metric blocks,
not one-electron Hessian approximations. Block-Jacobi PCG retains all coupling
in the full matrix-free action. Acceptance still requires the original true
residual, ||b-Gy|| <= sqrt(machine epsilon) ||b||. There is no tolerance change,
diagonal clipping, empirical shift, new iteration budget, or dense reduced
metric/inverse.

The inverse blocks use shared accepted-point storage O(sum_i d_i^2), not
O((sum_i d_i)^2), and are built only on the first nonzero inverse solve.
Full-AO analytic dual norms allocate no preconditioner blocks. No new AO-square
or integral copy is introduced per worker/rank. This is a constant-factor
preconditioner, not an asymptotic improvement or a universal spectral bound.

## Final matched MnF2 result

Job 280629, anode009; baseline/update sequentially in the same allocation,
same input, binding and numerical settings. GNU-time wall and process peak RSS
are from native runs, not profiling. BLAS uses one thread; nested teams disabled.

| Quantity | Baseline | Updated |
| --- | ---: | ---: |
| Total wall (s) | 39.21 | 24.95 |
| NEO metric wall (s), inclusive | 13.8586 | 1.2154 |
| HVP wall (s), inclusive | 21.0751 | 19.8776 |
| Accepted orbital iterations | 14 | 14 |
| HVP calls | 669 | 669 |
| Final total energy (Eh) | -1348.893353119308 | -1348.893353119308 |
| Peak RSS (MiB) | 514.18 | 511.91 |
| User + system CPU (s) | 963.19 | 735.88 |
| Process CPU / wall | 24.56 | 29.49 |

Total wall decreased 36.4%; the metric phase decreased 91.2%. RSS is effectively
unchanged: a two-MiB observation is not evidence of a systematic memory reduction.
The previous approximately 16-second metric observation was on another run;
the percentage above uses the matched 13.86-second baseline, not that historical
value. CPU/wall includes runtime spin/wait and is not a count of useful cores.

The final run performs 756 inverse solves and 26,064 metric actions, taking
1.2214 and 1.0592 seconds respectively (inclusive scopes overlap). Fifteen
accepted-point block constructions total only 0.0171 seconds. These counts
do not establish the old iteration count, which was not instrumented.

Preconditioner-only development pair, job 280621 on anode009:
wall 38.31 / 24.73 s, metric 13.7544 / 1.1648 s, both 14 steps and 669 HVP.
This isolates the main benefit before the HVP housekeeping changes.

## HVP audit: necessary actions versus redundant overhead

Final native low-level timers:

| Kernel | Calls | Wall (s) | Explanation |
| --- | ---: | ---: | --- |
| Exact AO-pair action | 688 | 6.5555 | 669 HVPs + 19 trial forwards |
| Fused AO-H1E forward/adjoint | 669 | 5.7741 | One fused sweep per HVP |
| AO-H1E adjoint | 727 | 4.1413 | 669 local responses + 43 structure adjoints + 15 gradients |

The later adjoint consumes the directional active-space gradient, which
depends on the earlier forward response. It is not the same RHS as the
adjoint already fused into the core response. Removing it as a duplicate
would remove a derivative term. The cached directional AO-pair image is also
already reused by the fixed-upstream two-electron pullback.

Two housekeeping changes are retained: CSR output rows are first-touched
and cleared by their owning workers, instead of a full serial clear; compact
stream accumulation and reduction share one team and reduce only actual
workers. Worksharing over CSR partitions also covers all rows when the runtime
supplies fewer workers than requested. Arithmetic within each CSR row remains
in the original edge order. Existing compact-stream private images are not
enlarged or claimed to have been eliminated.

These changes **do not demonstrate an MnF2 speedup**. The isolated job 280622
on anode009, preconditioner-only/final-housekeeping pair, gives total wall
24.56 / 24.79 s and HVP wall 19.5764 / 19.7067 s. Such small single-pair
differences do not distinguish an effect from timing variability;
the measured gain is from the mathematical metric preconditioner.

## Scaling and profiler qualification

Job 280624, anode017, same binary/input, native controls before sampling:

| Threads | Total wall (s) | HVP wall (s) | Metric wall (s) | CPU/wall | RSS (MiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 8 | 31.76 | 27.0319 | 1.1822 | 7.76 | 425.84 |
| 32 | 24.16 | 19.0140 | 1.2012 | 29.20 | 532.80 |

Both use 14 steps and 669 HVPs, with identical printed final energy. The
fourfold thread increase gives 1.31x total speed and 1.42x HVP speed; parallel
efficiency remains limited even though serial metric work is now small.
AO-pair/H1E contractions, rather than additional metric precision changes,
are the next evidenced target. The table used the eager-block development
binary; the final lazy change does not change the MnF2 inverse algorithm.

One gprofng capture preserves the exact development binary, experiment,
functions, call tree and quality report. Resolved kernels and barrier waits
agree qualitatively with the native timers. GNU 2.42 again warns that the
collection interval timer changed; percentages are not accepted as accurate
CPU shares, wall fractions, or estimates of removable barrier time.

## Validation and artifacts

Jobs 280621, 280622 and 280629 pass quotient/physical-map/Riesz tests, coupled
structure tests and shifted preconditioner tests at 1, 8 and 32 threads.
The final two jobs also pass exact active two-electron scalar/block and
CSR/compact tests. New analytic metric blocks match both the full coupled
action and independently implemented physical-feature norms for partial,
frozen and full supports. Certified inverse tests compare dense test-only
solutions over six RHS directions and scales 1e-8, 1 and 1e8, including full
OEO active-span, inactive-only and no-inactive cases.

Final parallel control array 280635, each baseline/update pair on the same
allocation, using the same NEO exact-integral settings:

| System | Steps, both | Final energy, both (Eh) | Wall baseline/update (s) | RSS baseline/update (MiB) |
| --- | ---: | ---: | ---: | ---: |
| 241 | 6 | -230.720590376550 | 1.89 / 1.85 | 940.12 / 956.83 |
| 240 | 8 | -343.446302348805 | 2.77 / 2.82 | 1680.56 / 1721.71 |
| FeCl2 | 3 | -2181.617644857086 | 11.91 / 11.74 | 436.26 / 430.31 |
| CERRAS | 6 | -397.093434948706 | 35.23 / 34.87 | 1988.00 / 2023.36 |

Small timing differences do not establish general speedups/regressions.
The control RSS is not uniformly lower; the largest observed increase is
41.15 MiB (2.45%) for 240. Additional inverse storage is orbital-local and
shared, not replicated by thread count; allocator/native RSS variation is
not equivalent to its live matrix payload. No heap attribution was collected
for that difference, so its cause is not asserted. FeCl2 records zero block
construction calls after the lazy change, unlike five unnecessary constructions
in the eager development version.

Local artifacts: results/studies/mnf2-metric-blocks-20261010/.
Remote: /home/guqqgroup/taoxia/xmvb-runs/mnf2-metric-blocks-20261010/.
Baseline, preconditioner-only and final binaries are preserved. The sampled binary
and its capture are under mnf2-riesz-check-280624/sampled. Final native hashes,
resources, CPU/node descriptions and phase logs are under mnf2-metric-280629.
The prematurely scheduled control array 280630 is excluded from comparisons.
