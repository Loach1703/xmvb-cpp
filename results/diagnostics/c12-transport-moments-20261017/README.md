# Fused finite-transport moments: C12 qualification

## Mathematical change

The finite-transport factor path previously made ten exterior generator
pullbacks per panel: four for each current-spin chart, and one for each origin
chart. Upper-chart moment pullback is linear. The four required contributions
for one current-spin chart are now accumulated in one occupation/shear walk
before a single pullback. Both orientations of the current-metric derivative,
both target factors, both origin terms and the selected transpose image remain.
Different spin topologies are not conflated.

After this fusion, independent off-diagonal upper-generator entries have
exclusive worker owners. The diagonal computes each row dot once and reduces
only orbital-length worker vectors. No carrier/factor panel is copied per
thread and no string-product/string-square storage is introduced. The factor
path retains the existing conservative16-column workspace budget; the new
simultaneous seed panels use at most7 such columns per spin. The single-term
and directional interfaces use the same common engine.

Theory: `article/transport_generator_moments.md`. The implementation changes
`exterior_transform.hpp/.cpp` and `transport_tiles.cpp`; independent tests are
in `test_representation_overlap.cpp`. Parent revision5f628433. The three
unrelated user-dirty files were not included or transferred.

## Native full-convergence results

Same original LOFLEA input, CAS(12,12), COV132, Ualpha=Ubeta924, HAO,
Libcint exact/cc-pVDZ, CIAH/Davidson. Hanhai22 anode40,32 physical cores,
affinity60-91, GCC14.1.0/CMake3.30.5, C++23 Release
`-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`, Eigen BLAS/LAPACKE and BLAS1.
Original thresholds and ITMAX2000; no optimization policy changes.

| Version / job | Steps | HVP | Final energy / Eh | SCF / s | Application / s | Peak RSS / MiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Saved scatter baseline /1264636 | 4 | 36 | -399.106646267331 | 28.207024 | 30.863387 | 2824.68 |
| Saved temporal physical fusion /1264713 | 4 | 36 | -399.106646267330 | 26.043647 | 29.009265 | 2725.68 |
| Saved shared-right moments /1264759 | 4 | 36 | -399.106646267333 | 25.845468 | 28.868012 | 2682.66 |
| Fused transport moments /1264773 | 4 | 36 | -399.106646267328 | 24.749263 | 27.623396 | 2679.75 |

Native application wall is4.3% lower than the immediate predecessor,4.8%
lower than1264713, and10.5% lower than the saved scatter baseline. These are
single-run measurements, not statistical speedup intervals. RSS is GNU-time
process maximum(KiB/1024), not Slurm's sampled accounting value. No major RSS
gain is claimed. The per-step HVP counts stay6,17,10,3; keyframes3,8,3,1.
The finite trajectories are not bitwise identical; final energies agree within
5e-12Eh. No old executable was rerun.

## Separate phase-instrumented runs

| Scope/counter | Saved1264759 | New1264773 |
| --- | ---: | ---: |
| Full application / s | 29.068800 | 26.449473 |
| Finite transport gradient / s | 2.685965514 | 0.593677933 |
| Exterior generator pullback calls | 440 | 176 |
| Exterior generator pullback wall / s | 2.317848341 | 0.190111607 |
| Exterior moment pullback wall / s | 0.004569729 | 0.002312342 |
| Complete coupled-image wall / s | 10.003442395 | 9.676410692 |

The exact10→4 calls/panel reduction explains440→176 calls. Generator wall
is12.2x smaller; the whole transport-gradient scope is4.5x smaller. Fusion
and output-owned parallelism were qualified together; no isolated attribution
to either one is claimed. Small triangular pullbacks were not the bottleneck.
All phase scopes overlap and are not additive time shares. Phase and native
run times must not be interchanged. No new profiler sampling claim is made.

## Correctness gates and provenance

Job1264773 completed with exit0 and9/9 regression targets passed. Tests include
selected nonclosed structure spaces, dense/minor cross-overlap references,
unequal alpha/beta spaces, one-column and7-column budgets,1/4-thread transport,
signed/zero/mixed-width fused moments, directional finite differences and
cooperative calls within an existing OpenMP team. Existing tolerances remain:
2e-12/3e-12 for action equivalence and2e-7 for the added central-difference
directional check. Full CIAH coupling and keyframe tests also passed. This
qualifies the changed contractions and this C12 convergence workload; it does
not assert universal speedups on untested molecules.

All four changed compiled source/test files match `sources.sha256`.
Application SHA256:
`a31106bd9fee7ce6f14cec8f794b669a2c8b12886af538fae2895a68981423dd`.
Input SHA256:
`092d545498678b0cf772c5c84b9eab674b694375558cd391a1b739b7bdcc35eb`.
Basis SHA256:
`d53e7e02d55fb99cdf823d50cf13fec375e0cbe1220ff61816aee7ab38e38ec9`.
The archive retains inputs, logs, traces, tests, resource records, flags,
source/binary hashes and node/affinity environment. Frozen binaries remain at
`/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015/c12-regression-1264773/bin/`
on Hanhai22 and are excluded from the repository archive.

`phase-summary.json` was generated with the existing
`tests/mcscf/summarize_c12_cost.py`. `comparison.csv` uses the existing saved-run
`row()` parser in `c12-temporal-responses-20261017/summarize.py` for the four
native rows above; it performs no calculation. The launcher is unchanged
`tests/mcscf/c12_regression.sbatch`.

## Next measured boundary

The physical coupled-image stage still costs9.68s over36 HVPs. The necessary
channel Gram/Hamiltonian products and selected-space output projections remain.
Future work must distinguish operations that can truly be shared from distinct
physical directions that require independent images; it must not discard CI
response or reinstate full-carrier/thread-private storage to gain speed.
