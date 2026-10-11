# Eight-edge CSR accumulation

## Scope

Single-node OpenMP/SIMD, exact AO integrals, NEO/Davidson. Baseline: 386ed52.
No integral replication, new scratch arrays, screening changes, relaxed solver
criteria, MPI, RI changes or molecule-dependent switches.

The saved gprofng capture from job 280820 resolves the runtime-selected AVX-512
`contract_rows` function as a hotspot. Its timer warning disallows quantitative
use of sample percentages. This round uses native phase times and matched runs
to assess the specific CSR accumulation change.

## Change

Each CSR output row already has one owner and shares its input coefficient
matrix with all workers. The former kernel accumulates four consecutive edges
before writing the output; the new kernel accumulates eight. Each scalar channel
still receives contributions in CSR edge order. Indices and integral values
are read once, not replayed for each channel tile. The final zero-to-seven edge
remainder is handled by the required scalar-edge SIMD loop.

This halves nominal output-row reads/writes in the main edge loop. It does not
halve total traffic: source-vector loads remain unchanged, as do CSR metadata
loads. Excluding metadata and assuming double precision, the source/output
traffic implied by the loop is 12 bytes per edge/channel before versus 10 after.
These are source-level access counts, not measured DRAM traffic. Compiler
register allocation and cache behavior determine the realized benefit.

Eight is an internal loop-unrolling choice, not a convergence parameter or
active-space threshold. The same source is runtime-selected for AVX-512, AVX2
and baseline ISAs. Native performance evidence below is from AVX-512 compute
nodes; no AVX2 speedup is claimed. The eight-edge loop uses `end - edge >= 8`
to avoid signed overflow near the 32-bit CSR edge-address limit.

This differs from the rejected experiment in job 280851: that implementation
revisited an entire CSR row for each eight-channel accumulator block. This
implementation never repeats the sparse metadata traversal.

## Verification

New tests cover long CSR rows, an empty row, edge remainders, overwritten
outputs and channel widths 0, 1, 7, 8, 9, 36 and 65, against a separately
assembled dense symmetric reference. The unchanged derivative tests also
compare scalar/block integral directions and fused HVP pullbacks.

Job 281230 passes six suites at 1, 8 and 32 threads (18 executions): NEO,
coupled response, spectral trust region, exact active integrals, AO-H1E
primitives and orthogonal direct CI. Job 281228 tests an initial version with
the same eight-edge arithmetic but the original additive loop-bound style;
job 281233 and array 281234 measure the final overflow-safe version.

All compilation and execution use Hanhai25 Slurm compute nodes. Native builds
retain `-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`, Eigen BLAS/LAPACKE, and
single-thread BLAS. Affinity is `OMP_PLACES=cores`, `OMP_PROC_BIND=close`,
`OMP_MAX_ACTIVE_LEVELS=1`; there is no wait-policy override. Inputs, executable
hashes, CPU placement, GNU time resources and phase logs are retained per run.

## Initial result

Job 281228, anode020. All runs converge in 14 steps with 669 HVPs. The maximum
matched final-energy difference is about 1.1e-11 hartree.

| Threads | Previous AO-pair (s) | Eight-edge AO-pair (s) | Previous wall (s) | Eight-edge wall (s) |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 9.048 | 8.498 | 26.704 | 26.098 |
| 16 | 5.012 | 4.944 | 17.316 | 17.216 |
| 32 | 2.986 | 2.820 | 13.729 | 13.700 |

The targeted improvement is modest and core-count dependent. Whole-program
changes are smaller, especially at 32 threads. Reduced source-level writes
alone are not proof of a substantial end-to-end speedup or better scaling.

## Final matched runs

Job 281233, anode033; matched compact CPU sets at each thread count. All runs
converge in 14 steps with 669 HVPs. The maximum matched energy difference is
about 1.6e-11 hartree. Timings are single native captures, not repeated medians.

| Threads | Previous AO-pair (s) | Eight-edge AO-pair (s) | Previous wall (s) | Eight-edge wall (s) | RSS, previous / new (MiB) | New average CPU cores |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 7.822 | 7.518 | 23.489 | 23.088 | 430.9 / 449.1 | 7.60 |
| 16 | 4.254 | 4.007 | 15.144 | 14.936 | 437.4 / 449.8 | 14.71 |
| 32 | 2.615 | 2.652 | 12.717 | 12.774 | 498.6 / 493.5 | 28.46 |

The final target improves by 3.9% and 5.8% at 8 and 16 threads, but regresses
by 1.4% at 32. There is **no stable 32-thread MnF2 speedup** across the two
allocations: the initial target gain was 5.6%. Whole-program gains are at most
a few percent and peak RSS does not consistently decrease. No additional
allocation has been introduced; observed RSS variation is not evidence of an
algorithmic memory reduction. CPU accounting includes spin waits, not only
useful work. Relative-to-eight-thread 32-thread efficiency is 45.2% before
and 45.2% after; single-thread efficiency was not measured.

Array 281234 runs five independent matched 32-thread allocations. Compare
within each pair, not against job 281233. All exit statuses are zero; steps and
orbital HVP counts are unchanged. The maximum matched energy difference is
5e-12 hartree.

| System | Steps | HVPs | AO-pair, previous / new (s) | Wall, previous / new (s) | RSS, previous / new (MiB) | New energy (hartree) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 241 | 6 | 10 | 0.129 / 0.115 | 1.989 / 2.036 | 922.1 / 934.5 | -230.720590376551 |
| 240 | 8 | 13 | 0.337 / 0.339 | 4.091 / 4.021 | 1722.5 / 1769.5 | -343.446302348805 |
| FeCl2 | 3 | 302 | 1.080 / 0.968 | 7.732 / 7.637 | 428.2 / 407.9 | -2181.617644857090 |
| MnF2 | 14 | 669 | 2.793 / 2.683 | 13.325 / 13.023 | 507.6 / 496.2 | -1348.893353119306 |
| CERRAS | 6 | 30 | 1.240 / 1.142 | 39.469 / 39.404 | 1999.6 / 1973.9 | -397.093434948705 |

The target phase improves by 10.7%, 10.3%, 3.9% and 7.9% for 241, FeCl2,
MnF2 and CERRAS; 240 is effectively unchanged (0.4% slower). The small kernel
reordering is retained on this limited evidence, but is not a major parallel
improvement. Particularly, a faster target phase does not imply a faster total:
241 total wall time grows by 2.3%. No molecule-specific selector is added to
hide a less favorable measurement.

Further work should focus on shared integral traversals and eliminating
carrier/transform copies, not additional unrolling constants. The underlying
ERI arithmetic and source traffic remain unchanged; this patch cannot remove
those costs or establish near-ideal thread scaling.

Raw evidence lives under `results/studies/locality-20261011/` in
`mnf2-storage-281228/`, `neo-verify-281230/`, `mnf2-storage-281233/` and
`simd-run-281234/`. `locality_results.py` parses the native MnF2 logs.
