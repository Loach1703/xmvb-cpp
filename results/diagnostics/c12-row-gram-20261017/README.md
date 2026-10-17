# Row-partitioned channel Gram: C12 qualification

## Change

The necessary oriented Gram is G=A^T B, with exact exchange weights already
in B. Partitioning its common row index gives G=sum_t A_t^T B_t. Unlike an
output grid where each small product scans all rows, every worker reads its
own row segment for a full local product. The consumed input halo supplies
disjoint partial-output storage. A barrier precedes the exclusive-entry
reduction and another precedes scratch reuse. Neither live channel slab is
overwritten or copied. Independent orbital and CI cotangents remain separate.

The available team-size partial Grams must fit the existing halo. Otherwise
workers directly own output blocks of the same Gram, with no partial array
or reduction. This prevents a tiny halo from creating scalar-product/barrier
explosion; it does not select a different CI representation or wave function.
No new explicit numerical workspace or larger64MiB cap is introduced. BLAS
packing remains library-owned as before. See `article/channel_gram_rows.md`.

Parent revision b159c9e5. Changed production file: `src/mcscf/ci/sigma_tiles.cpp`.
Independent tests: `tests/vbscf/unit/test_orthogonal_direct_ci.cpp`. Three
unrelated user-dirty files were neither staged nor transferred.

## Same-workload native convergence

LOFLEA, CAS(12,12), COV132, Ualpha=Ubeta924, HAO, exact Libcint/cc-pVDZ,
CIAH/Davidson from the original input guess. Hanhai22 anode40,32 physical
cores, affinity60-91, GCC14.1.0/CMake3.30.5, C++23 Release
`-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`, Eigen BLAS/LAPACKE, OpenMP32,
BLAS1. Original thresholds and ITMAX2000 retained; no strategy changes.

| Version / job | Steps | HVP | Final energy / Eh | SCF / s | Application / s | Peak RSS / MiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Saved scatter baseline /1264636 | 4 | 36 | -399.106646267331 | 28.207024 | 30.863387 | 2824.68 |
| Saved transport-moment predecessor /1264773 | 4 | 36 | -399.106646267328 | 24.749263 | 27.623396 | 2679.75 |
| Initial row-Gram candidate /1264784 | 4 | 36 | -399.106646267331 | 21.325744 | 24.243956 | 2704.70 |
| Final row-Gram geometry /1264786 | 4 | 36 | -399.106646267332 | 21.394834 | 24.205419 | 2650.95 |

The final native application is12.4% faster than its saved immediate predecessor
and21.6% faster than the saved scatter baseline. The two new candidates are
not statistical replicates:1264786 added the focused full-square test and
replaced narrow-halo scalar partial tiles by direct output ownership. Both
share C12's full-row Gram geometry. No old binary was rebuilt or rerun.
RSS is GNU-time process maximum (KiB/1024); no large memory reduction is
claimed. The final value is2.59GiB, versus2.62GiB in the predecessor. HVP
counts6,17,10,3 and keyframes3,8,3,1 are unchanged, as are36 physical CI
responses. Native and instrumented final energies agree within a few1e-12Eh.

## Separate phase-instrumented attribution

| Scope / counter | Saved1264773 | Initial1264784 | Final1264786 |
| --- | ---: | ---: | ---: |
| Full application / s | 26.449473 | 24.308098 | 24.682620 |
| Complete coupled images / s | 9.676410692 | 7.743862966 | 7.989972917 |
| Channel compute / s | 5.540055210 | 3.212042977 | 3.235936756 |
| Gram scope / s | 2.661782951 | 0.551000969 | 0.563850898 |
| Gram GEMM products | 50688 | 25344 | 25344 |
| Generated sweeps | 129 | 129 | 129 |
| Generated batches | 1171 | 1171 | 1171 |
| Source elements read | 1120445436 | 1120445436 | 1120445436 |

The final Gram scope is4.7x faster and complete coupled-image scope17.4%
shorter. These are overlapping scopes, not additive time shares or profiler
sample percentages. Instrumented runs are separate from native timings above.
The source reader counter is unchanged: the reduction concerns internal Gram
operand rereads, not fewer external callback reads. Logical input reuse is
derived from the matrix decomposition; no DRAM-counter measurement is claimed.

Each C12 Gram uses32 full78x78 partials. Their1,557,504 bytes alias the dead
input halo and must not be added to the64MiB workspace. The recorded live
workspace payload is57,277,200 bytes for these sweeps, unchanged. Different
summation order is expected; physical state space, accuracy and outputs are
unchanged.

## Gates and evidence

Job1264786 completed with exit0,9/9 regression targets passed. Added focused
complete8-orbital3e/3e tests use width2 and2/4 workers, distinct response/state
carriers, two temporal requests and parities0,+1,-1. They compare every image
and both integral cotangents to independent non-tiled `apply` and
`integral_adjoint` operators at the existing3e-11 tolerance. The fixture's
6272-double halo holds four36-square partials. Existing tests cover unequal
spins, nonclosed support, zero state/directions, directional and forcing-only/
adjoint-only consumers, nested teams, and now explicit four-worker minimal
halos including full/vacuum spin strings. Full coupled HVP/keyframe tests pass.
All changed compiled production/test hashes match the final local files.

`evidence.tar.gz` retains both new jobs' raw input, output, CMake flags,
convergence trace, tests, time/RSS, node/affinity and source/binary hashes,
excluding frozen binaries. Binaries remain on Hanhai22 under each job's
`/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015/c12-regression-JOB/bin/`.
`phase-summary.json` was generated with `tests/mcscf/summarize_c12_cost.py`;
`comparison.csv` uses the existing saved-run `row()` parser in
`c12-temporal-responses-20261017/summarize.py`. Neither parser runs workloads.
The unchanged launcher is `tests/mcscf/c12_regression.sbatch`.

The remaining coupled action still costs7.99s, including necessary channel
Hamiltonian products and selected-space output/generator contractions. The
next target is to identify shared output projections or directional products
without dropping the nonzero-Ritz-residual normalization/gauge terms.
