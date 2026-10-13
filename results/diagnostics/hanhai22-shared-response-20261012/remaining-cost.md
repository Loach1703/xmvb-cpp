# Remaining cost audit — 2026-10-13

## Scope and evidence

Production source is unchanged from b9d775e4. No threshold, coordinate, HVP,
or convergence-rule changes. No old executable was rerun. All workload execution
was through Slurm on Hanhai22, 32 OpenMP threads, single-threaded BLAS, core
binding, frozen executable from build 1250646. Input/binary hashes, flags,
environment, node and CPU are recorded in each run directory.

| Job | Workload | Outcome |
|---|---|---|
| 1250649 | Existing QOCGUM scalar trace, reused | Converged, 4 macros, 52 coupled expansions |
| 1250663_0 | New LOFLEA full scalar trace | Converged, 4 macros, 34 coupled expansions |
| 1250663_1 | LOFLEA software VTune attempt | Collection denied by ptrace policy; no valid profile |
| 1250665_1 | LOFLEA two-macro phase trace | Expected application exit 2 at diagnostic iteration limit |

The last run is attribution-only, not a convergence test. VTune failure is
retained in `remaining-cost-1250663/LOFLEA-sampled/run.err`; no machine security
settings were changed. Results below are instrumented scope timings, not a new
production speed comparison or sampled function percentages. Existing native
LOFLEA wall time remains 74.030271 s (`keyframe-lean-comparison.csv`).

## Full-run scalar scopes

Automatically extracted with `summarize-remaining-cost.py`; raw tables in
`remaining-cost.csv`. These scopes are not a complete additive time budget.
Coupled expansion includes action and response-space observation; it is not
the exclusive orbital-HVP kernel time.

| Scope | QOCGUM / s | LOFLEA / s |
|---|---:|---:|
| Structure eigensolves, initial plus accepted points | 3.741 | 47.751 |
| Coupled expansion action | 2.512 | 21.252 |
| Expansion metric images | 0.074 | 1.865 |
| Analytic diagonal construction | 1.514 | 2.842 |
| Internal keyframes | 1.147 | 3.816 |
| Preconditioner applications | 0.024 | 0.004 |
| Small AH Ritz solves | 0.007 | 0.002 |

Further optimization of small AH diagonalizations is not the priority here.

## Important new localization: structure-action transforms

LOFLEA has 12 active electrons in 12 active orbitals, 132 VB structures, and
924 unique strings for each spin. The orthogonal carrier is a 924 by 924
spin-string product tensor. Davidson remains in the selected VB-structure
space; this is not a new determinant-space eigensolver.

In the two-macro run, three structure eigensolves cost 34.954257 s. All direct-CI
sigma calls across the entire run cost 12.078576 s (368 calls). Hence at least
22.875681 s of eigensolve time lies outside sigma. This is a lower bound, not
an exclusive timing of the exterior transforms: structure expansion/projection,
transforms, Davidson bookkeeping and other overhead must still be separated.

Code facts:

- `src/vbscf/structures/assembly/action.cpp:1198`: each direct-CI structure action
  expands the structure block, applies two forward transforms, applies sigma,
  applies four adjoint transforms to the H and overlap images, then contracts
  back to structures.
- `src/vbscf/structures/orthogonal_ci/exterior_transform.cpp:178`: transforms
  use serial shear/pair loops. Left transforms update rows of a column-major
  matrix, resulting in strided memory access. Right transforms have contiguous
  column updates but are also serial.
- `src/vbscf/structures/orthogonal_ci/sigma.cpp`: sigma itself already has a
  bounded-workspace OpenMP implementation. Its process CPU/wall ratio is about
  30 in this trace, versus about 13 for the whole trial-forward scope. These
  ratios are utilization clues, not proof of useful arithmetic occupancy.

Proposed first implementation, after adding exclusive transform timers:

1. Left transforms: assign disjoint columns or column tiles to workers; retain
   the full ordered shear sequence within each tile. Improve column locality.
2. Right transforms: assign disjoint row ranges; apply ordered shears to the
   corresponding contiguous column segments.
3. Preserve reversed shear order for adjoints. Do not parallelize dependent
   shears. Avoid one barrier per shear and avoid full-carrier thread copies.
4. Measure forward, adjoint, structure expansion/contraction and sigma separately.
   Verify transform/adjoint identities and full structure actions, then compare
   converged new runs against existing native data. No speedup is claimed yet.

## Confirmed smaller redundancies

### Retain converged selected H/S images

`src/core/eigensolver.cpp:480` already constructs both H and overlap images of
the converged roots. The return path retains overlap images but not H images.
`src/vbscf/derivatives/hessian/context/response_cache.cpp:106` recomputes both
when selected-state images are absent. Passing both through the accepted-point
context avoids one repeat structure action per used point. Storage is only
structure count times selected-state count, not a spin-product-sized cache.
Retain the actual H image; substituting eigenvalue times overlap image would
discard the finite Davidson residual.

### Do not build an unused terminal preconditioner

`src/vbscf/optimization/backends/neo.cpp:693` builds the next chart and its
analytic diagonal before checking final convergence. Four accepted steps
therefore cause five diagonal constructions. Final constructions cost about
0.30 s in QOCGUM and 0.55 s in LOFLEA and do not lead to another step.
Separate geometry/step mapping from curvature construction; avoid building a
second chart or repeating SVD work. This is useful but not the main LOFLEA cost.

### Reuse metric images returned by coupled actions

The structure coupling action returns a metric image that the backend discards,
while expansion obtains the same-vector metric image separately. A combined
image interface could share this result without retaining stale cross-point
data. This is secondary to structure transforms; do not reintroduce full norm
inversions to achieve it.

## Lower-priority algorithmic option

Davidson retains its basis and H/S images within each fixed orbital point, but
only the selected roots survive across macro points. A bounded structure-space
correction subspace could be recycled across points. New-point H/S images must
be recomputed; old images are invalid. Evaluate total RHS work, not merely
block-call counts, because warm-start blocks also cost actions. This is a
hypothesis, not a demonstrated saving, and should follow transform optimization.

## What this audit does not justify removing

- H/S structure action itself: required for the selected-space eigenproblem.
- Coupled orbital–structure response: the phase trace has zero nested
  `structure_solve` calls; this run is not repeatedly solving a full outer
  response equation inside every HVP.
- Complete-space physical residual checks: the Davidson callback is gated
  near convergence, not called for every expansion. LOFLEA's selected covalent
  space does not enable the complete-active-space callback.
- Correctness fixes or the existing outer convergence thresholds.
