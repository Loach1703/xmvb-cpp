# CIAH preconditioner comparison

Array job 1247579, Hanhai22, 32 physical cores per task, exact libcint
integrals, Davidson, identical frozen executable SHA-256
`9a0a7196af4264ff3f09d326be86cb98cbc4568473da586bd9cf2d0f85304863`.
Implementation commit: `a8d60c8c`. No recompilation or optimizer change.

Five array tasks run concurrently. Each task evaluates Hessian-diagonal and
then one-electron on the same node/allocation and input. The outer thresholds
are gradient 1e-3 and energy change 1e-7 Eh. The script does not override the
input/default iteration limit. The Slurm allocation has a one-hour resource
limit, not an extra 30-step optimizer limit. BLAS uses one thread;
OpenMP uses 32 bound cores, with nested parallelism disabled.

These labels refer to the current CIAH **signed** analytic/one-electron
curvature preconditioners. This is not a comparison of the historically
different projected-NEO solver, nor of block-LBFGS's positive initial inverse.
The physical orbital step map and exact Hessian action are the same in both
variants.

## Completed pairs

| System | Preconditioner | Steps | Orbital HVP | Final energy (Eh) | Wall (s) | Peak RSS (MiB) |
|---|---|---:|---:|---:|---:|---:|
| F2 | Hessian-diagonal | 3 | 24 | -198.751155830525 | 0.31 | 17.19 |
| F2 | One-electron | 3 | 23 | -198.751155830527 | 0.27 | 16.96 |
| 241 | Hessian-diagonal | 3 | 26 | -230.720590392853 | 2.67 | 802.66 |
| 241 | One-electron | 3 | 23 | -230.720590392890 | 2.56 | 830.57 |
| 240 | Hessian-diagonal | 3 | 40 | -343.446302353194 | 5.49 | 1563.07 |
| 240 | One-electron | 3 | 34 | -343.446302353166 | 5.33 | 1534.58 |

Final energies agree within 3.7e-11 Eh for these pairs. One-electron uses
fewer HVPs in all three. The small native timing differences are single-run
observations, not statistically established speedups; initialization is
included in native wall time.

## Incomplete pairs at this snapshot

| System | Completed Hessian-diagonal: steps / HVP / wall / RSS | Energy (Eh) | One-electron |
|---|---|---:|---|
| MnF2 | 11 / 169 / 5.56 s / 551.84 MiB | -1348.893353223710 | Still running, over 200 steps |
| FeCl2 | 20 / 429 / 11.91 s / 385.86 MiB | -2181.624464318884 | Still running, over 170 steps |

These two one-electron runs were not cancelled and have no final timing/RSS
or converged energy yet. Do not turn missing results into zeros or speedups.
FeCl2 remains near the previously observed high-energy branch; MnF2 has
slow energy/gradient progress. Both final endpoints remain to be checked.

Thus analytic Hessian-diagonal is **not uniformly superior**. It has a large
practical advantage in the two transition-metal runs so far, while
one-electron has fewer actions in the three small molecular runs. This pilot
does not justify a molecule-specific dispatch rule or a universal ranking.

`summarize_ciah_preconditioners.py RUN_DIRECTORY` parses output, native
GNU-time records and NEO counters into CSV automatically. Incomplete action
totals count only recorded, accepted iterations. Raw artifacts are retained
under `ciah-preconditioners-1247579/` on the cluster and in the local results
mirror; node/CPU/affinity, input/binary hashes and exact commands are saved.
