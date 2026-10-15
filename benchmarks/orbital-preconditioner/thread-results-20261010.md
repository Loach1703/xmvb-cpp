# Complete 241 NEO thread scaling

Hanhai25 job 280481, anode006. Cartesian cc-pVDZ, exact libcint, HAO,
CAS(6e,6o), FULL structures, Davidson, NEO, one-electron preconditioner.
Identical input and production-optimized phase-instrumented binary from
`xmvb-phase-audit-20261010/updated` for all native runs. Independent Fe-porphine
convergence jobs were not stopped. No numerical algorithm was modified.

All three native runs converged in six iterations to -230.72059037655 Eh;
differences are approximately 1e-12 Eh. Gradient/energy tolerances are 1e-3
and 1e-7 Eh. One full run per thread count, in order 1, 8, 32, in one allocation.
This is a shared node, not an exclusive reservation: initial system load
averages were 216.61, 222.90, 232.32 on a 256-core node. First-run cold-cache
effects and other workloads have not been removed; these are exploratory
single-run measurements, not statistical speedup claims.

| Threads | GNU wall (s) | User + system CPU (s) | Peak RSS (MiB) | Speedup vs 1 thread |
|---:|---:|---:|---:|---:|
| 1 | 12.22 | 12.18 | 812.7 | 1.00 |
| 8 | 2.34 | 15.45 | 770.8 | 5.22 |
| 32 | 2.06 | 46.85 | 954.2 | 5.93 |

## Native phase evidence

Internal wall clocks exclude a little process startup/shutdown overhead.
Trial forward and gradient are subsets of SCF, not additive to SCF. Chart
contains its projection and factorization scopes; nested times must not be
summed as independent stages.

| Phase (s) | 1 thread | 8 threads | 32 threads |
|---|---:|---:|---:|
| AO integral provider | 2.770356 | 0.540015 | 0.314770 |
| AO integral input | 0.414320 | 0.170831 | 0.145197 |
| Full input preparation | 3.192415 | 0.715198 | 0.463888 |
| Full SCF iterations | 8.950570 | 1.575041 | 1.531969 |
| Seven trial forward evaluations | 1.944284 | 0.316683 | 0.412518 |
| Seven trial gradient completions | 0.919873 | 0.133967 | 0.084476 |
| Chart (15 calls) | 0.055991 | 0.056713 | 0.057347 |

The 8-to-32-thread SCF improvement is only 2.7%; most of the end-to-end
improvement comes from integral preparation. Trial forward becomes slower
in this particular pair while gradient completion improves. The chart is
essentially thread-independent but only about 57 ms at 32 threads, so it does
not explain the bulk of the 1.53-second SCF time. Total CPU consumption grows
threefold from 8 to 32 threads with little SCF wall improvement. This supports
investigating parallel overhead/saturation; it does not alone distinguish
load imbalance, spinning, bandwidth, and serial regions.

## Sampling evidence and limitations

VTune software sampling was denied on anode006 by the kernel ptrace scope;
the native runs succeeded, but job 280481 ended with status 1 on the separate
profiler step. No security settings were changed.

GNU gprofng job 280483 completed, with the same six-iteration convergence and
energy. The capture lasted 2.174 seconds. Its resolved hot symbols include
`gomp_barrier_wait_end`, `gomp_team_barrier_wait_end`,
`apply_exact_ao_pair_kernel`, `apply_ao_h1e_fused`, `apply_ao_h1e_transpose`, and
`apply_ao_h1e`. This provides useful priorities: OpenMP synchronization and
exact AO-pair/J-K response contractions, not the tiny orbital block factor.

However, gprofng reports `Collection interval timer period was changed
(10007 -> 0); profile data may be unreliable`. Therefore the displayed
sampling percentages are not certified runtime fractions, and barrier samples
are not proof of equivalent removable wall time or a specific load imbalance.
Do not claim that the program could be accelerated by the sum of barrier
percentages. A same-node wait-policy ablation and narrower kernel wall timers
would distinguish spinning from actual useful work and isolate the limiting
parallel phase before changing numerical code.

Raw native logs and reports are under
`results/studies/241-thread-scaling-20261010/`; full remote captures are retained
under `/home/guqqgroup/taoxia/xmvb-runs/241-thread-scaling-20261010/`.
