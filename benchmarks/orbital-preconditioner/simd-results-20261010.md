# Shell workspace and direct-CI layout optimization

## Changes and ownership

The immutable libcint evaluator, normalization and optimizers are now shared
within an OpenMP process. Each worker owns only reusable shell outputs and
libcint scratch. Exact quartet evaluation, raw-quartet prescreening, RI metric
and RI three-center evaluation use the modern caller-cache interfaces. The
previous returning-by-value two-electron APIs and per-worker evaluator copies
are removed, not retained as alternate production paths. One-electron shell
output allocation is unchanged in this round.

`LibcintShellWorkspace` owns scratch and three shell buffers, not complete AO
integral arrays. Returned references must be consumed before their buffer is
reused. Smoke comparisons explicitly copy snapshots when they need two blocks
at the same time. There is no shared mutable scratch or thread-local singleton.

Direct-CI intermediates store string/RHS entries contiguously within each
orbital-pair column. This transposes storage, not the physical contraction:
sigma and integral-adjoint GEMMs are algebraically transposed together.
The beta gather is a contiguous SIMD loop. Sigma reduction buffers are first
allocated by their owning workers and reduced over the actual team size.
Adjoint buffers use the actual number of available strings rather than always
allocating 160 rows. Coefficient packing and the exact selected-structure
projection are unchanged. No additional determinant-layer calculation is added.

This is a single-node OpenMP/SIMD change, **not MPI support**. Full private AO
matrices in one-electron scatter and compact AO-pair action still require a
separate output-ownership design. A single irreducible tile can also exceed an
existing workspace budget. Do not infer that all large-space memory risks are
resolved.

## Execution and numerical checks

Hanhai25, partition CPU-256C768GB, one MPI task, 32 OpenMP workers,
single-threaded BLAS, pinned cores. Release flags remain
`-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`, Eigen BLAS enabled. No host-specific
AVX/AVX-512 flags or fast-math are added. Baseline production binary is from
the preceding exact-block revision; input and binary hashes are saved.

Build/test jobs 280532, 280540, 280542 and 280551 passed. Checks at 1, 8 and
32 threads include:

- Independent legacy libcint quartet reference over all 81 s/p/d shell
  combinations, including zero-valued components and Cartesian normalization.
- Shared optimizer calls with worker-owned scratch; repeated equal-sized
  quartets preserve output and cache pointers.
- RI normalization, positive/singular metric rank and factor contractions.
- Existing direct-CI explicit matrix references, arbitrary selected-structure
  projection, distinct spin spaces, integral adjoints and response pullbacks.
- Added 210-string, two-RHS diagonal-spectrum reference spanning a partial
  final tile, action/adjoint agreement, dynamic and nested OpenMP teams.

Job 280552 additionally ran converged F2 RI calculations: both variants use
five iterations and yield -198.750976531946 hartree. GNU wall time is 0.16 s
for both; RSS is 71.8/64.0 MiB. This small case establishes end-to-end RI
correctness, not a meaningful timing improvement. The printed AO provider time
is zero on that RI path and cannot be used to attribute RI construction cost.

## Native end-to-end results

Final layout: array 280543; each independent case has baseline then updated
in one allocation. NEO, exact libcint, Davidson, projected gradient tolerance
1e-3 and energy-change tolerance 1e-7 hartree. No short iteration/time cutoff.
GNU process wall time includes startup/output; RSS is process peak in MiB.
Shared-node single paired runs do not establish statistical significance.

| Case | Steps, both | Wall baseline/new (s) | RSS baseline/new (MiB) | Mean occupied cores baseline/new |
| --- | ---: | ---: | ---: | ---: |
| 241 | 6 | 1.86 / 1.78 | 966.9 / 944.1 | 23.0 / 23.8 |
| 240 | 8 | 3.01 / 3.03 | 1747.7 / 1680.1 | 19.1 / 19.1 |
| FeCl2 | 3 | 12.22 / 12.28 | 437.3 / 440.3 | 25.4 / 25.6 |
| MnF2 | 14 | 133.33 / 132.67 | 553.2 / 523.3 | 8.7 / 8.8 |
| CERRAS | 6 | 39.61 / 38.67 | 2030.4 / 1972.7 | 25.3 / 25.2 |

Mean occupied cores are GNU user-plus-system CPU seconds divided by GNU wall
seconds, including possible spinning; they are not useful-work core counts.
Final energies agree within 3e-12 hartree in these cases. AO provider times
baseline/new are 0.298253/0.270417 (241), 0.483637/0.449951 (240),
0.278497/0.227787 (FeCl2), 0.441548/0.223206 (MnF2),
0.601261/0.556488 seconds (CERRAS).
Overall improvement is modest; faster integral construction does not imply
faster SCF. Neither 240 nor FeCl2 has a measured total-time gain here.
MnF2's roughly nine occupied cores out of 32 allocated indicate poor overall
parallel resource use, not a specific proven cause. Its long wall time is
already present in the baseline; this change does not resolve that bottleneck.

The first candidate (280535) exposed small-adjoint leading-dimension overhead;
it is retained as evidence, not used as the final-layout timing table.

## Isolated contraction evidence

Job 280544 uses an identical public-API benchmark harness linked to each
variant's library on one compute node. Cases have 20 or 462 strings per spin
and block widths 1 or 4. Sigma is warmed once then timed for five calls;
adjoint is timed for one call. Setup and full VBSCF costs are excluded.

For 462 strings at 32 workers:

| Kernel | Block width | Baseline per call (ms) | Updated per call (ms) |
| --- | ---: | ---: | ---: |
| Sigma | 1 | 9.704 | 6.111 |
| Sigma | 4 | 26.603 | 14.895 |
| Integral adjoint | 1 | 9.436 | 5.233 |
| Integral adjoint | 4 | 12.375 | 10.915 |

Checksums agree to floating-point precision. These measurements do not
establish a universal improvement: at one worker the 462-string, one-RHS
adjoint is slower (46.305 to 52.537 ms); at eight workers it also regresses in
this capture. Tiny kernels are sensitive to team startup and timing noise.
The four-RHS large kernel improves at all three measured thread counts.

GCC's saved vectorization report confirms the new beta-gather loop uses
16-byte vectors. The evidence supports contiguous SIMD and BLAS layout, not
an AVX-512 claim. Full reports and unlinked diagnostic objects are retained.

## Allocation evidence

Heaptrack job 280549, full 241 calculation, eight workers. Both profiled
variants converge in six iterations to the same printed energy.

| Metric | Baseline | Updated |
| --- | ---: | ---: |
| Allocation calls | 2,653,901 | 599,899 |
| Temporary allocations | 804,103 | 244,850 |
| Reported peak heap | 1.01 GB | 1.00 GB |

Allocation calls fall about 77.4%, but live-heap peak changes little.
The baseline report resolves over one million allocation calls from
`CINT2e_drv`, with quartet evaluator call paths. Caller-owned scratch removes
that repeated allocation mechanism. Profiled runtime/RSS include collector
overhead and are not native speed or RSS evidence. Reported heap and RSS may
differ because allocation extent, residency and tool accounting differ.

## Next performance target

Use these changes as foundations, not a completed deep parallel overhaul.
The next measured target remains exact AO-pair and AO one-electron response
scatter: assign output tiles, bound private reduction storage, and establish
thread scaling before adding MPI. MPI needs explicit integral/data ownership
and communication planning; simply multiplying ranks would replicate caches.

## Reference implementation sources

- [libcint caller-cache integral driver](https://github.com/sunqm/libcint/blob/master/src/cint2e.c)
- [PySCF string scatter/GEMM/gather](https://github.com/pyscf/pyscf/blob/master/pyscf/lib/mcscf/fci_contract.c)
- [PySCF shell-block direct J/K](https://github.com/pyscf/pyscf/blob/master/pyscf/lib/vhf/nr_direct.c)

These are implementation references, not benchmark baselines. No qcint GPL
SIMD source is copied into the project.

Raw logs/captures are under results/studies/simd-shell-20261010 locally and
/home/guqqgroup/taoxia/xmvb-runs/simd-shell-20261010 remotely. Benchmark scripts
and the optional `test_orthogonal_direct_ci --benchmark` harness reproduce the
experiments without running computations on local or login nodes.
