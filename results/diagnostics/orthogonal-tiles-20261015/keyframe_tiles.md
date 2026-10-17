# Shared fixed-state and finite-transport contractions

## Implementation and ownership

- `structures/assembly/states.*` owns fixed-state normalization, the small
  projected state Hamiltonian, the complete overlap multiplier and the
  original-frame coefficient gradient. The gradient evaluator no longer
  duplicates these operations.
- `prepare_active_adjoint` optionally emits its H/M images from the same sigma
  sweep. Keyframes first obtain M C, normalize the prescribed state span,
  then collect H C and unshifted adjoints together. The complete small
  Hamiltonian supplies the overlap multiplier afterward. Off-diagonal state
  terms are retained; no state is assumed stationary at a keyframe.
- `structures/assembly/transport_tiles.cpp` owns the finite-transport VJP.
  Four formerly independent chart contractions and the reference CI image
  share factor generation and seeds. Both current-metric derivatives and the
  reference-direction connection remain present.
- Factor panels obey the existing 64 MiB overlap allowance. A conservative
  sixteen-alpha-plus-beta-panel bound controls width. Original touched
  support is sampled directly, with no selected-space closure assumption.
  Workers own rows of shared panels, not replicas. Complete broad support
  retains a five-image moment order only when it fits the same allowance.
- The main orbital and independent structure direction sweeps are **not**
  fused in this round: their separate images still feed the coarse
  preconditioner. That larger interface change is deliberately left explicit.

There are no changes to CIAH decisions, physical gradient/HVP formulas,
preconditioner, spin space, convergence thresholds, metric-solve criteria,
or integral screening. No new rank/amplitude approximation or certification.
These changes remove redundant work; they do not change the exponential
complete-CAS carrier scaling.

## Native COV full convergence

Original LOFLEA HAO guess, exact libcint cc-pVDZ, CAS(12,12), singlet, 132 COV
structures. NEO/CIAH, structure Davidson, analytic Hessian-diagonal,
gradient 1e-3 and energy change 1e-7 Eh. All runs below are on anode50,
eight actual OpenMP cores in one LLC, BLAS one thread. Old binaries were not
executed again. Process time is GNU time end-to-end, not a single SCF step.

| Version / job | Steps | HVP / response | Final energy (Eh) | SCF (s) | Process (s) | RSS (GiB) |
|---|---:|---:|---:|---:|---:|---:|
| Saved local selected-metric / 1263333 | 4 | 36 / 36 | -399.106646267332 | 67.355929 | 71.48 | 2.510 |
| Previous transport fusion / 1263716 | 4 | 36 / 36 | -399.106646267331 | 60.480704 | 65.18 | 2.506 |
| Shared fixed-state images / 1263892 | 4 | 36 / 36 | -399.106646267332 | 52.317433 | 57.08 | 2.522 |
| Plus finite-transport VJP / 1263912 | 4 | 36 / 36 | -399.106646267331 | 50.668144 | 55.19 | 2.505 |

Final process wall is 15.33% below the preceding round, 22.79% below the saved
local 71.48 s run. That saved local run is **not** the GitHub baseline.
Macro microstep counts (10,31,17,3), HVP counts (6,17,10,3), keyframes
(3,8,3,1), and rejected trials (1,1,0,0) are unchanged. Input SHA-256 is
5d3086bf489db0d039c33197dc3ca61f3e5dc85bc4e1a1aaada292383ef0a858.

These are single shared-node samples, not averages. The previous native run
used cores184--191; both new native variants used136--143. Unmodified stages
also vary in measured time, so the entire wall-time difference is not assigned
to one eliminated function. Exact work-count changes below identify the
removed operations independently of machine variation.

## Scoped phase and work evidence

Saved phase run1263717 versus new1263913, same anode50. Nested timers are not
additive and are not sampled CPU percentages.

| Operation | Before | After |
|---|---:|---:|
| Ordinary CI sigma sweeps | 88 | 77 |
| CI response sweeps | 88 | 88 |
| Exterior-map constructions | 432 | 234 |
| Ordinary overlap adjoints | 37 | 4 |
| Directional overlap adjoints | 11 | 0 |
| Independent overlap actions | 26 | 15 |
| New finite-transport sweeps | 0 | 11 |
| Keyframe transport gradient (s) | 4.153584628 | 2.384617430 |
| Complete keyframes (s) | 12.912470661 | 9.322409675 |

There are more product projections (3240 ->3570): collecting the adjoint's
already-generated image requires a projection, but removes a separate sigma
and its expansion. The optimization does not claim every counter decreases.
In the existing `ciah_keyframe_active_cost` record, `action` now includes the
joint action/adjoint sweep and `adjoint=0` means no separate adjoint replay;
the physical derivative has not been omitted.

## FULL single-step check

| Version / job | HVP | CI actions | After-step energy (Eh) | Process (s) | RSS (GiB) |
|---|---:|---:|---:|---:|---:|
| Saved GitHub94693457 /1262623 | 8 | 39 | -399.322429054264 | 93.77 | 3.373 |
| Saved later local /1263111 | 8 | 27 | -399.322366416760 | 80.84 | 3.370 |
| New /1263914 | 8 | 27 | -399.322366416770 | 72.34 | 3.383 |

These FULL tasks have226512 structures and deliberately stop at one accepted
step (`max_iterations`, application exit2). They are **not converged** runs.
The new trajectory agrees with the later local trajectory at1e-11 Eh; the
different remote trajectory predates the earlier CI/coordinate changes and
is not attributed to this round. Native FULL time includes its large final
structure-weight output, unlike the shorter application timer. RSS remains
similar, with a small increase relative to the saved local FULL sample.

## Verification and artifacts

Final build/test job1263904 passed8/8 tests: representation overlap, CI
transport, selected CI transport, keyframe, structure tiles, sparse orbital
quotient, F2 HAO/OEO Hessian diagonal. Expanded tests cover fused versus
separate fixed-state images/adjoints, full state multipliers, independent
occupied-minor derivatives, arbitrary nonclosed selections, equal/unequal
spins, empty beta space, dense/sparse support, one-column/wide budgets and
one/four-thread scatters. Existing keyframe orbital/coefficient finite
differences and transported HVP tests remain active.

The intermediate fixed-state-only build1263885 also passed8/8. Native and
phase data are generated with the existing Python parsers, not manually
transcribed. `keyframe_tiles_native.csv`, `keyframe_tiles_phases.csv`, JSON,
source hashes, accounting and the evidence archive retain provenance. The
theory note `article/ci_transport_tiles.tex` was expanded and its four-page
PDF compiled and inspected. Unrelated local RI changes remain untouched.
