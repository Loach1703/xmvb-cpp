# MnF2 phase scaling: 8, 16, and 32 cores

## Measurement

Revision `5ace1f4`, exact Libcint, NEO/Davidson, unchanged gradient and
energy tolerances (1e-3 and 1e-7 hartree). Slurm job `281244` ran each
configuration once on `anode020`, in one 32-core allocation. These are
single-run diagnostic timings, not statistically established speedups.
Compilation, regression tests, and profiling all ran on compute nodes.
No production algorithm or tolerance was changed for this measurement.

OpenMP uses compact affinity; BLAS uses one thread. The 8/16/32-core sets
span respectively 1/1/2 NUMA domains and 1/2/4 L3 domains. Therefore the
32-core comparison also changes memory locality; hardware memory-bandwidth
counters were unavailable. The allocation was not node-exclusive.

| Threads | Application wall, s | Steps | Orbital HVPs | Energy, hartree | Peak RSS, MiB | Average CPU cores | Efficiency relative to 8 cores |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 26.689 | 14 | 669 | -1348.893353119326 | 448.2 | 7.63 | 100.0% |
| 16 | 17.017 | 14 | 669 | -1348.893353119312 | 473.6 | 14.75 | 78.4% |
| 32 | 13.598 | 14 | 669 | -1348.893353119318 | 491.1 | 28.26 | 49.1% |

Average CPU cores are GNU-time user plus system CPU divided by elapsed
time. They include OpenMP spinning and are not useful-work core counts.
Peak RSS increases by 42.9 MiB from 8 to 32 cores; no large-array copying
was introduced in this diagnostic round.

## Phase results

| Phase, seconds | 8 cores | 16 cores | 32 cores | Speedup, 16 to 32 |
| --- | ---: | ---: | ---: | ---: |
| AO-pair action | 9.117 | 4.822 | 2.808 | 1.72 |
| Fused AO one-electron response | 7.271 | 3.915 | 2.380 | 1.65 |
| AO one-electron transpose | 4.634 | 2.520 | 1.532 | 1.65 |
| Three AO phases combined | 21.022 | 11.257 | 6.720 | 1.68 |
| Everything outside these three AO phases | 5.666 | 5.760 | 6.879 | 0.84 |
| NEO projected solve | 2.815 | 2.812 | 3.208 | 0.88 |
| Projection factor/model solve | 1.322 | 1.349 | 1.442 | 0.94 |
| Projection residual/certification | 1.493 | 1.463 | 1.766 | 0.83 |
| Physical metric solve | 1.296 | 1.295 | 1.374 | 0.94 |
| Selected-state directional action | 0.107 | 0.123 | 0.297 | 0.41 |
| Active backward | 0.075 | 0.092 | 0.219 | 0.42 |
| Structure solve | 0.041 | 0.073 | 0.321 | 0.23 |

Only the three AO scopes are combined here: they do not overlap one
another. Other rows are nested diagnostics and must not be added together.
In particular, projection factor/residual are inside projected solve;
physical metric solve also contributes to residual certification.

Doubling cores saves 3.419 s overall, rather than the ideal 8.509 s.
Of the 5.090 s gap to ideal halving, 1.091 s comes from the three AO
phases and 3.998 s from everything outside them. This is accounting, not
proof that the latter work is all serial. The projected solve alone is
23.6% of the 32-core application wall time.

## Load balance and profiler qualification

| Worker-loop occupancy indicator | 8 cores | 16 cores | 32 cores |
| --- | ---: | ---: | ---: |
| AO-pair | 99.2% | 98.5% | 95.8% |
| Fused AO one-electron | 97.6% | 97.3% | 95.0% |
| AO one-electron transpose | 96.0% | 96.0% | 96.8% |

The indicator compares accumulated worker loop durations against the
slowest-worker capacity of each invocation. It measures loop imbalance,
not useful CPU utilization or exact barrier wait. These results do not
support a large AO load-balance tail as the dominant remaining problem.
Balanced loops can still be limited by memory bandwidth, NUMA access,
cache pressure, or common runtime overhead; this capture does not separate
those causes.

Companion job `281246` profiled the same executable at 32 cores on the
same node. It converged in 14 steps to -1348.8933531193 hartree. Both perf
hardware and software-event probes were denied. GNU gprofng reported
`Collection interval timer period was changed (10007 -> 0)`; its percentages
are therefore not used as accurate time fractions. Its function names
identify AO kernels and OpenMP barrier routines as investigation leads,
consistent with, but not a quantitative replacement for, native phase timers.

## Source-level interpretation and next work

1. `optimization/neo/response_solver.cpp` reconstructs whitening, the
   projected structure pseudoinverse, and the relaxed-model spectrum at
   each subspace iteration. The measured factor phase does not scale with
   OpenMP threads. Audit dependency changes before reusing any factors;
   a previous version already attempted structure-inverse caching, so
   repeating that change without a measured avoided-work count is not useful.
2. Residual certification calls the physical dual norm, whose iterative
   metric inverse is in `optimization/trust_region/retraction.cpp`.
   Its current strict true-residual certification must remain unchanged.
   Identify repeated right-hand sides, projections, or metric images before
   proposing batching or reuse; do not substitute a weaker certificate.
3. Selected-state, active-backward, and structure-solve stages regress at
   32 cores. The direct-CI code in `structures/orthogonal_ci/sigma.cpp`
   creates teams for small kernels. Measure invocation sizes and team
   counts before fusing regions or assigning bounded worker counts.
   Use shared outputs and bounded scratch, not per-thread large replicas.
4. The AO kernels still matter, but the next large total-time opportunity
   is outside them. Additional SIMD unrolling alone cannot remove the
   measured non-scaling floor. Cross-NUMA affinity is a separate experiment,
   not an explanation established by these data.

## Reproduction and evidence

Run `mnf2-balance.sbatch SOURCE_DIR` on Hanhai25. The companion capture is
`mnf2-clock.sbatch SOURCE_DIR 32`; the optional argument selects the capture
thread count within the existing allocation.

Parse completed logs with `locality_results.py` and `phase_results.py`.
Local evidence: `results/studies/locality-20261011/mnf2-balance-281244/`,
`mnf2-clock-281246/`, `scaling-281244.csv`, and `phases-281244.csv`.
Remote evidence is under
`/home/guqqgroup/taoxia/xmvb-runs/native-20261011/`.

Executable SHA256:
`1c2fb68686f7333158aba37e21a12a2cf3243f353cb003cc1ce8ff40098d4f2b`.
Input SHA256:
`79034c0fdb5baee2432e81d90c7316de2b0bf2edfe35166b184f2c28f761a822`.
The sparse exact-two-electron and RI AO-H1e HVP regression executables both
passed in the native compute-node job.
