# AO-H1E pair-group contraction

## What transfers from PySCF

The [AO2MO driver](https://github.com/pyscf/pyscf/blob/master/pyscf/lib/ao2mo/nr_ao2mo.c)
assigns independent output rows to an OpenMP team and uses small reusable buffers
for two-stage BLAS transforms. The [direct J/K driver](https://github.com/pyscf/pyscf/blob/master/pyscf/lib/vhf/nr_direct.c)
partitions shell blocks, shares reordered density inputs and consumes each
integral block for multiple operators. Its private JK buffers are bounded but
can still be large per thread; that policy is not appropriate to copy blindly.
These are execution/data-layout techniques, not orthogonal-orbital optimization
assumptions. Sources inspected on 2026-10-11.

Our active-pair CSR contraction already owns disjoint output rows, shares its
inputs and uses runtime-selected SIMD kernels. Replacing that screened sparse
graph with a dense AO integral tensor just to call GEMM would discard screening
and increase memory. The useful transfer for this round is **reuse within an
integral block**, applied to the AO-H1E response rather than another small CI
parallel region.

## Reordering the AO-H1E contraction

Unique ERI references are already ordered by CSR storage. Contiguous references
with the same left AO pair therefore reuse their left AO indices, column offsets,
diagonal symmetry factor and Coulomb source entry. The scalar forward, transpose
and fused actions now process such groups inside the existing resource-derived
dynamic chunks.

For a fixed left pair, its Coulomb contributions to that same output position
are summed in one scalar and written once at the end of the group. Exchange and
right-pair contributions are still accumulated immediately. Coincident indices
do not invalidate the separation: those terms are additions to the same output,
not dependencies of the Coulomb sum. This changes floating-point sum order but
not the operator or its transpose. The original symmetry factors are preserved.

One compile-time forward/transpose kernel replaces the former scalar per-integral
accumulators. The independent multi-RHS kernel remains required and unchanged.
Compact streams may be unordered: grouping uses actual consecutive left-pair
equality, not an assumption that rows are sorted. No alternate numerical fallback
or decoded-stencil cache is added.

The arithmetic complexity remains linear in the number of stored ERIs for each
scalar response. This is removal of repeated index decoding, source loads and
Coulomb output updates, not a new asymptotic integral algorithm. New scratch is
only scalar local state; the existing AO-sized worker accumulators remain. No
thread-private integral or active-pair image is introduced.

## Verification and measurement

All builds and executions use Hanhai25 Slurm compute nodes, exact Libcint
integrals, NEO/Davidson, the original inputs and original convergence criteria.
Release flags are `-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`, Eigen BLAS/LAPACKE
enabled, single-thread BLAS, `OMP_PLACES=cores`, `OMP_PROC_BIND=close`, and
`OMP_MAX_ACTIVE_LEVELS=1`; no wait-policy override or precision relaxation.
The immutable baseline is the preceding revision eb28dd8.

The existing gprofng capture (job 280820) identified these response loops as
hotspots. Its timer warning prevents quantitative use of sampling percentages;
the comparisons here use native phase timers and GNU time resources.

Regression job 281214 completes all six suites at 1, 8 and 32 threads (18
executions). The exact AO-H1E tests now include two and five AO basis functions,
mixed pair orientations, unordered compact streams, coincident indices, empty
integral lists, workspace reuse, dynamic teams and nested calls. Scalar and fused
results are compared to an independently assembled dense operator and its
transpose at the existing absolute tolerance of 2e-13. Existing block response,
active integral, direct-CI, coupled-response and NEO tests also pass.

An initial native experiment, job 281213, tested index grouping and Coulomb sums.
The final implementation also explicitly holds the fixed Coulomb source entries
in scalars; job 281215 measures that final code. Executable/input hashes and
affinity records are retained for each run. Job array 281217 performs matched
32-thread checks on five molecules in parallel.

### MnF2 thread comparison

Job 281215, anode033; each matched pair uses the same compact allocated CPU
set. Every run converges in 14 steps with 669 orbital HVPs. Final energies span
-1348.893353119316 to -1348.893353119312 hartree.

| Threads | Previous wall (s) | Grouped wall (s) | Previous RSS (MiB) | Grouped RSS (MiB) | Grouped average CPU cores |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 24.749 | 23.723 | 412.4 | 430.8 | 7.60 |
| 16 | 16.054 | 14.873 | 456.3 | 484.1 | 14.74 |
| 32 | 13.481 | 13.097 | 501.5 | 504.6 | 28.30 |

| Threads | Previous fused / transpose (s) | Grouped fused / transpose (s) | Combined response reduction |
| ---: | ---: | ---: | ---: |
| 8 | 7.622 / 4.730 | 6.626 / 4.188 | 12.4% |
| 16 | 4.147 / 2.621 | 3.511 / 2.230 | 15.2% |
| 32 | 2.612 / 1.627 | 2.305 / 1.458 | 11.2% |

The target improves at all three core counts. End-to-end reductions are 4.1%,
7.4% and 2.9%, respectively. This is not an RSS reduction: observed peaks are
slightly higher, although the implementation adds no dynamically allocated
scratch or physical array. Native RSS variability does not identify its cause.
Average process CPU includes OpenMP spinning; it is not useful-work accounting.
Eight-to-32-thread efficiency is 45.3% for the grouped code versus 45.9% for the
baseline, despite faster absolute times. No single-thread efficiency is claimed.

### Five molecular checks

Job array 281217 uses five independent allocations, one matched baseline and
grouped run per case. These CPU placements differ from job 281215, so timings
must not be compared between jobs. All exits are zero, steps/HVPs are unchanged,
and the largest matched final-energy difference is 3e-12 hartree.

| System | Steps | HVPs | Wall, previous / grouped (s) | SCF, previous / grouped (s) | RSS, previous / grouped (MiB) | Grouped energy (hartree) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 241 | 6 | 10 | 2.121 / 2.026 | 1.017 / 0.931 | 920.5 / 937.4 | -230.720590376550 |
| 240 | 8 | 13 | 4.145 / 4.256 | 1.523 / 1.493 | 1784.6 / 1679.4 | -343.446302348805 |
| FeCl2 | 3 | 302 | 7.884 / 7.756 | 7.281 / 7.150 | 427.8 / 408.1 | -2181.617644857093 |
| MnF2 | 14 | 669 | 13.708 / 13.048 | 13.078 / 12.401 | 480.4 / 516.2 | -1348.893353119315 |
| CERRAS | 6 | 30 | 39.554 / 39.487 | 36.766 / 36.633 | 2039.1 / 1988.0 | -397.093434948706 |

The combined fused/transpose phase decreases by 11.8%, 9.5%, 13.0%, 13.2% and
14.8%, respectively. The SCF stage decreases for every case; small total-time
differences need cautious interpretation. In particular, 240 total time grows
because input preparation increases by about 0.14 s, while its SCF stage falls
by 0.03 s. Input preparation is not modified by this patch. CERRAS receives
little end-to-end benefit because AO-H1E is a small fraction of its total cost.

Raw logs, native resources, input/binary hashes and affinity records are saved
under `results/studies/locality-20261011/mnf2-storage-281213/`,
`mnf2-storage-281215/`, `neo-verify-281214/` and `simd-run-281217/`.
`locality_results.py` produces `ao-group-initial.csv` and `ao-group-final.csv`.
All jobs complete successfully; no calculations ran locally or on login nodes.

## Dependency boundary

The early AO-H1E forward and transpose responses already share an ERI pass.
The late transpose consumes adjoints produced by structure/active response and
cannot be moved into that early pass. This patch preserves that dependency.

A later shared sweep with the independent active-pair directional contraction
would need to prepare its input before the early response and preserve
disjoint active-pair output ownership. A naive unique-integral scatter would
instead replicate the full active-pair image per worker, defeating the memory
goal. No such fusion is claimed here.
