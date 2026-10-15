# Exact response and OpenMP workspace audit

## Scope

All builds and calculations ran through Slurm on Hanhai25. Job 280509
built and tested the updated source; array 280511 ran the three independent
cases concurrently, with baseline and updated binaries sequentially in each
allocation. Each calculation used 32 OpenMP threads, single-threaded BLAS,
exact integrals, Davidson, NEO, gradient tolerance 1e-3 and energy tolerance
1e-7 hartree. Inputs and binary hashes are recorded in each run directory.

## Changes

- Integral-only directional batches reuse accepted half-transformed integrals
  without computing unused directional half-transforms.
- Resident directional images are applied in memory-budgeted blocks; streamed
  integral-only batches share accepted row tiles. Streamed fixed-adjoint
  requests still use separate directional sweeps.
- Move, rather than copy, the accepted half-transform into its cache.
- AO one-electron forward/transpose workspaces are initialized by their owning
  workers. Output reduction stays in the same OpenMP region and uses the actual
  team size, preserving ascending thread summation order.
- Batch one-electron reduction writes directly into its returned matrices,
  removing two full temporary sum matrices.

The initial batch-only experiment (280490) did not show a clear end-to-end
speedup: the main NEO coupling path uses the scalar API. The second experiment
includes one-electron initialization/reduction changes on that main path.

## Completed matched results

Times below are GNU time total process wall time, including startup/output;
RSS is GNU time peak resident memory in MiB. These are single paired runs,
not statistically established speedups.

| Case | Steps, both | Final energy, both (hartree) | Baseline wall (s) | Updated wall (s) | Baseline RSS (MiB) | Updated RSS (MiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 241 | 6 | -230.720590376551 | 1.96 | 1.91 | 959.3 | 949.6 |
| 240 | 8 | -343.446302348805 | 3.05 | 3.02 | 1831.8 | 1813.2 |
| FeCl2 | 3 | -2181.617644857089 | 11.93 | 11.65 | 445.1 | 436.4 |

Program-internal SCF times, excluding startup/output, were respectively
1.460934/1.387690, 2.163406/2.157274 and 11.580967/11.274458 seconds
(baseline/updated). Printed final energies and iteration counts agree.
Small RSS differences should not be interpreted as proven elimination of
the dominant memory allocations.

Job 280509 passed active two-electron directional/block regressions and AO
one-electron forward/transpose/fused/block regressions at 1, 8 and 32 threads,
as well as NEO globalization. Coverage includes CSR/compact integrals,
resident/streamed representations, dynamic and nested OpenMP teams.

## Wait-policy ablation

Job 280503 compared the batch-only binary with default OpenMP waiting and
OMP_WAIT_POLICY=PASSIVE/GOMP_SPINCOUNT=0. Program wall times were:

| Case | Default (s) | Passive (s) |
| --- | ---: | ---: |
| 241 | 1.944282 | 1.964071 |
| 240 | 3.047576 | 4.733246 |
| FeCl2 | 11.664118 | 14.802051 |

Passive waiting is not adopted. Reduced spinning CPU time is not equivalent
to reduced wall time. Earlier gprofng timer warnings also preclude treating
its barrier percentages as exact elapsed-time fractions.

## Remaining memory and parallelism constraints

1. OpenMP workers share the immutable AO integral graph; it is not copied per
   worker. However, AO one-electron scatter still uses full thread-private AO
   matrices, two per worker in the fused action. First-touch and parallel
   reduction do not remove that replication.
2. Compact AO-pair integral scatter similarly retains thread-private output
   images. CSR row-owned output avoids that particular replication.
3. The 64 MiB directional block budget accounts for compact private images
   and temporary transforms, but is not an absolute process memory cap:
   returned images and an irreducible single direction can exceed it.
4. Next changes should establish output/tile ownership and bounded scatter
   reduction, while measuring the cost of any repeated integral scans. Do
   not trade reduced storage for unmeasured repeated work.
5. SIMD should target contiguous active-pair/RHS loops; it does not require
   copying the integral graph. Layout and write ownership precede forced
   vectorization of scatter loops.
6. No MPI implementation is added here. Running multiple unmodified ranks
   would replicate the integral graph and caches. A useful MPI design needs
   explicit integral ownership, node-shared read-only data where appropriate,
   bounded local tiles and communication of vector/tile results rather than
   complete replicated operators.

Raw evidence: results/studies/exact-block-20261010/block-run-280511,
block-run-280490 and wait/wait-policy-280503 (local, not committed).
