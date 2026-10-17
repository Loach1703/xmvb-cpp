# C12 projected-channel experiment: no production path retained

All six numerical candidates passed 11/11 regressions and converged in
4 accepted steps / 36 orbital HVPs / 36 structure responses. None beat the
saved native control. **All experimental production changes were removed.**
The numerical sources again match 17a4825d (whose retained implementation is
d54c9095). There is no projected-image API, extra dispatch, candidate toggle,
or numerical fallback left in the maintained implementation.

Retained work is the contraction derivation in `article/channel_support.md`,
a bounded arithmetic-ownership diagnostic in `audit_channel_support.cpp`,
and this evidence. Unrelated local RI/ignore edits were neither transferred
nor staged. Rejected source is only the noncompiled evidence patch.

Restored-state job **1265047** rebuilt all relevant targets and passed **11/11
tests in 3.51 s**. It ran no native SCF control. The saved control above was
not rerun. No jobs remain active at the final accounting query.

## Matched end-to-end convergence

Hanhai22, anode42, CPU affinity 96-127, 32 physical cores, GCC14.1/C++23
production Release flags. OpenMP close/cores; BLAS threads=1. The molecular
input is **LOFLEA HAO, CAS(12,12), 132 covalent structures, cc-pVDZ, exact
Libcint, 150 basis functions, ISCF=9/Davidson**. It is not a separate carbon-C12
geometry. Original guesses, precision, response and optimization policies are
unchanged. The old native control was read from saved data, not rerun.

| Variant / job | Native application / s | SCF / s | Peak RSS / MiB |
|---|---:|---:|---:|
| Saved retained / 1264965 | 17.764194 | 15.287966 | 2720.47 |
| Projected boundary / 1265025 | 19.368058 | 16.662374 | 2704.68 |
| Row dot and channel-index reuse / 1265032 | 18.266775 | 15.568900 | 2712.52 |
| Fused Gram and boundary images / 1265037 | 18.094699 | 15.425786 | 2705.23 |
| Also project ordinary H/S / 1265038 | 18.450621 | 15.861753 | 2695.23 |
| Remove redundant boundary writes / 1265039 | 18.212159 | 15.528575 | 2650.33 |
| Cyclic boundary ownership / 1265046 | 19.164215 | 16.502251 | 2746.91 |

All native per-step HVP counts are 6,17,10,3; final energies range from
-399.106646267332 to -399.106646267330 Eh. CSV retains unrounded values.
Timing differences are single-run observations, not statistical confidence
intervals. RSS is GNU-time process peak, not attributed heap; its variation
does not establish a memory reduction. Every candidate kept one shared
64 MiB channel slab lease and used only small extra row/metadata workspaces.

The failed first build 1265023 was a test-local identifier collision, fixed
before any native measurement; it is not a numerical failure or timing run.

## What was verified mathematically

The explicit projected-image contract permits K[J,J] on halo-only rows,
without assuming closure of a selected VB structure space. Full sigma still
requires K[:,J]. Derivative consumers were restricted only after auditing
their source support and upper-generator pullbacks. The complete covalent
boundary certificate emits exact channels from deficient-prefix intervals.

Tests cover arbitrary selected supports, unequal alpha/beta populations,
temporal requests, mixed complete/projected image requests, integral adjoints,
and independent gradient/HVP and triangular-pullback checks. Ordinary H/S
projection was separately tried and removed. No physical response was dropped.

Dense interior products and sparse boundary products share tile storage.
The fused version consumes the boundary Gram before overwriting its source,
reuses its emitted indices for base/directional images, and reuses the existing
small Gram reduction scratch. These identities are exact, but the resulting
small indexed products plus partition/synchronization did not yield a native
speedup at this workload. This experiment does not show that every possible
implementation of the identities must be slower.

## Arithmetic balance is not measured parallel efficiency

Census 1265040 streams bounded 4096/32768-row diagnostic tiles. The native
production tile partition is not claimed identical to these diagnostic tiles.
For COV the boundary work model is sum(k_r^2)=23083410:

| Diagnostic rows | Continuous 32-worker capacity | Cyclic capacity |
|---|---:|---:|
| 4096 | 51082656 | 27751392 |
| 32768 | 33908384 | 24623168 |

Capacity is 32 times each tile's maximum assigned arithmetic work, summed
over tiles. This predicts better arithmetic balance under cyclic ownership,
not cache behavior, wall time, or measured CPU efficiency. The actual cyclic
candidate was slower and was rejected. The census retains at most one tile
of widths plus fixed-size statistics, not a global product-space table.

## Profile and evidence

1265033 captures the frozen 1265032 candidate using gprofng on compute nodes.
Quality reports no collection errors/warnings. Resolved leads include the
E-transpose gather callback, channel scatter, GEMM, AO contractions and OpenMP
waits. CPU sampling shares are not wall-time percentages; waiting samples do
not alone prove load imbalance. No hardware memory/cache counters were captured.

Native timing/RSS, phase diagnostics and profiled data are separate. In fused
candidates the recorded Gram scope includes boundary image computation, so it
cannot be interpreted as a directly comparable pure-Gram timer. Final cyclic
phase logs explicitly count `gram_fused_boundary_rows`.

`evidence.tar.gz` contains native/phase logs, traces, resource records,
hashes/toolchains, census and profiler reports, and the rejected source patch;
executables and the large raw sampling capture are excluded. Frozen binaries
remain in remote `c12-regression-JOB/bin/`. The raw sampling experiment remains
at `/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015/c12-cost-1265033/loflea-profile/capture.er`.
No new speedup is claimed for this round.
