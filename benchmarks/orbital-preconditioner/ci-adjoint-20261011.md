# Shared direct-CI adjoint tiles

## Implementation

The integral adjoint previously accumulated a full active-pair-square output
per worker, followed by a serial reduction. The replacement scatters the left
and right spin-string operands into two shared tiles and assigns disjoint
two-dimensional blocks of the output to workers. The contraction and its
derivative are unchanged; only accumulation order changes.

Tiles cover RHS blocks, alpha strings and beta strings. Their combined explicit
operand payload is bounded by 64 MiB, independently of the thread count. The
required shared output and BLAS packing are **not** included in this limit.
The output grid uses the most nearly square exact factorization of the requested
worker count. OpenMP work sharing covers all blocks even when a dynamic or
nested team contains fewer workers. Scatter and contraction barriers prevent
input overwrite while another worker is consuming a tile.

Two independent one-body generator requests in the directional structure
action now share one parallel region. Inputs are views, not concatenated carrier
copies; the two outputs were already required. The remaining dependent generator
action stays separate. On MnF2, generator phase entries fall from 2007 to 1338
for the same 669 HVPs. This does not establish a timing improvement: their
measured cumulative time was 0.138 s before and 0.141 s after.

No MPI, replicated full integrals, numerical tolerance changes, or old adjoint
fallback were introduced. The required shared dense active-pair output remains;
this change does not remove all fourth-power active-integral storage.

## Independent kernel measurement

Hanhai25 Slurm job 280880; one compute-node allocation; compact CPU placement;
single-thread BLAS; matched Release builds (`-O3 -DNDEBUG`, Eigen BLAS/LAPACKE,
OpenMP). The synthetic space has 48 orbitals, 48 alpha strings and 48 beta
strings, **one electron of each spin**, and 1176 packed orbital pairs. It is not
a 48-electron molecular calculation. Time covers one integral-adjoint call;
RSS covers the entire benchmark process. Each run verifies its adjoint against
the action bilinear form, with contracted value 2.0003329281006925.

| Threads | Previous adjoint (s) | Shared adjoint (s) | Previous RSS (MiB) | Shared RSS (MiB) |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 0.05418 | 0.04124 | 154.8 | 127.3 |
| 16 | 0.06967 | 0.04858 | 239.4 | 128.8 |
| 32 | 0.11261 | 0.06969 | 410.8 | 146.8 |

At 32 threads the measured call is 38% shorter and process peak RSS is 64%
lower. However, both implementations become slower beyond eight threads in
this small test. This is an output-memory improvement, not evidence of ideal
thread scaling. Single-call timings are exploratory, not statistical estimates.

## Matched molecular regression

All molecular calculations ran on Hanhai25 compute nodes with exact integrals,
unchanged inputs and convergence criteria. Baseline numerical code is revision
9943444 (the immutable baseline executable contains equivalent numerical code
plus phase instrumentation). Job 280881 used identical compact placements for
each matched MnF2 pair. All variants converged in 14 steps with 669 HVPs.

| MnF2 threads | Previous wall (s) | Shared wall (s) | Previous RSS (MiB) | Shared RSS (MiB) | Shared average CPU cores |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 26.922 | 27.486 | 462.6 | 438.5 | 7.66 |
| 16 | 17.376 | 17.485 | 466.1 | 475.5 | 14.85 |
| 32 | 14.673 | 14.638 | 497.9 | 500.3 | 28.68 |

The final 32-thread energy is -1348.893353119312 hartree. This round establishes
no material MnF2 end-to-end speedup. Eight-to-32-thread speedup is 1.88, or 47%
of ideal scaling relative to eight threads. Average CPU occupancy at 32 threads
is 90%; occupancy alone does not imply useful parallel efficiency. No matched
single-thread run was made in this round.

Job array 280882 ran five independent cases in parallel, each with a matched
32-thread baseline and updated run. Placements are not identical to job 280881;
compare only within each matched pair.

| System | Steps, both | Previous wall (s) | Shared wall (s) | Previous RSS (MiB) | Shared RSS (MiB) | Updated energy (hartree) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 241 | 6 | 2.181 | 2.238 | 965.1 | 932.8 | -230.720590376550 |
| 240 | 8 | 3.566 | 3.590 | 1692.4 | 1728.1 | -343.446302348805 |
| FeCl2 | 3 | 7.117 | 7.004 | 408.7 | 397.5 | -2181.617644857099 |
| MnF2 | 14 | 15.003 | 14.777 | 510.4 | 513.4 | -1348.893353119326 |
| CERRAS | 6 | 36.855 | 37.403 | 1991.1 | 2012.3 | -397.093434948704 |

All exits are zero. The maximum matched final-energy difference is approximately
1.1e-11 hartree. End-to-end changes are small and do not demonstrate a speedup.
Peak RSS does not decrease for every molecule: the bounded shared operand tiles
can be larger than the previous small operands, even though replicated outputs
have been removed.

Job 280879 passed six regression suites at 1, 8 and 32 threads (18 executions).
Coverage includes exterior derivatives, integral-adjoint finite differences,
directional gradient/HVP checks, tile boundaries, dynamic teams and nested teams.

Raw logs: `results/studies/locality-20261011/ci-adjoint-280880/`,
`mnf2-storage-280881/`, `simd-run-280882/`, and `neo-verify-280879/` under the
same parent. `locality_results.py` produces `ci-final.csv` from job 280881 logs.
The dedicated kernel driver and Slurm script are `ci-adjoint.cpp` and
`ci-adjoint.sbatch`; the old reference source is a benchmark artifact only.

## Remaining AO work

The three measured 32-thread AO phases still total about 7.72 s:
active-pair directional action 3.07 s, fused AO one-electron response 2.79 s,
and late transpose response 1.86 s. This round does not modify those kernels.

The early one-electron forward and transpose responses are already fused.
Sharing their integral traversal with the independent active-pair directional
transform requires reorganizing preparation and consumers. The late transpose
response depends on newly computed structure/active adjoints; it cannot simply
be merged with the early traversal. That dependency must be preserved in any
future shared-sweep optimization.
