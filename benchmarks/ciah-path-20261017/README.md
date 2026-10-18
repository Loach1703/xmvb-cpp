# Cached advancement of clipped CIAH proposals

The historical results below describe the original path-advancement change.
The2026-10-18 driver additionally admits a remaining path only after its
cached images satisfy the existing forcing condition for the current RHS.
This corrects the C12 FULL five-step regression without a threshold change;
see `../c12-c14-channels-20261018/README.md` and
`article/ciah_clipped_path.md` for its derivation and native comparisons.

## Change

An admitted AH proposal can require several component-limited pieces. Retain
that direction and its existing H image while following those pieces; do
not re-solve the AH equation or expand Davidson merely to take another piece.
Any kept true keyframe invalidates the unused path. After completion or
invalidation, reproject the current complete RHS before further expansion.

`CiahStepPath` borrows the existing `CiahResult` images and stores only a
pointer and a remaining fraction. It does not copy physical arrays, create
a new Hessian representation, or allocate per-thread workspaces. Clipped or
replayed pieces are not reported as converged AH roots for the current RHS.
The physical macro acceptance, coordinate/CI transport, HVP, preconditioner,
all thresholds, and microstep/keyframe budgets are unchanged.

Derivation: `article/ciah_clipped_path.md`. This is distinct from the failed
positive-interior and split-image experiments in
`article/ciah_subspace_steps.md`; those branches remain removed.

## Verification

Build/test job **291713** on Hanhai25 passed all seven focused tests:
`ciah`, `ci_transport`, `selected_ci_transport`, `keyframe`,
`coupled_structure`, `hessian_diagonal`, `orbital_step_geometry`.
The CIAH test includes exact segmented-image reconstruction with a nonidentity
metric, unchanged HVP count, unclipped identity, keyframe invalidation,
current-RHS reprojection and subsequent outside-space residual enrichment.

Native array **291714_[0-4]** ran all five inputs to convergence, without an
added SCF step limit. Exact Libcint integrals, CIAH/Davidson, analytic Hessian
diagonal, gradient threshold `1e-3`, energy-change threshold `1e-7 Eh`;
32 physical cores, cores/close binding, BLAS threads 1. Compilers/modules:
GCC 14.3.0 and CMake 3.31.9; C++23, Release `-O3 -DNDEBUG -g1
-fno-omit-frame-pointer`. Existing qualified Libcint/OpenBLAS libraries were
reused. No numerical task or build ran locally or on a login node.

| Input | Steps | HVP | Final energy / Eh | GNU wall / s | Peak RSS / MiB | Replayed pieces |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 240 | 3 | 27 | -343.446302545194 | 2.33 | 1743.07 | 1 |
| 241 | 4 | 29 | -230.720590533669 | 1.47 | 996.82 | 17 |
| FeCl2 | 13 | 170 | -2181.624464871226 | 3.68 | 444.67 | 86 |
| F2 | 3 | 18 | -198.751155830547 | 0.35 | 60.00 | 9 |
| MnF2 | 11 | 160 | -1348.893353291661 | 2.66 | 506.48 | 22 |

The final binary SHA256 is
`0bdd0eae85ae78b4c2b0696329063a9a55b9e06ca8d57c5a8578b3eada6c761a`.
Raw outputs, source/input/binary hashes, affinity, environment, per-macro and
per-piece traces, exit status and GNU time are in
`results/diagnostics/ciah-path-20261017/borrow/`.
`summary.csv` was generated from these logs and the saved reference logs,
not from manual trajectory transcription. `source` names the base commit
plus this patch; unrelated pre-existing RI worktree changes were excluded
from the execution mirror and left untouched.

## Existing references, not rerun

| Input | Saved version / job | Cluster / cores | Steps | HVP | GNU wall / s | Peak RSS / MiB |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| 241 | c491230b / 291576 | Hanhai25 / 32 | 4 | 44 | 2.03 | 984.30 |
| 240 | 48c145ba / 1264434_0 | Hanhai22 / 32 | 4 | 30 | 4.53 | 1569.29 |
| FeCl2 | 48c145ba / 1264383 | Hanhai22 / 8 | 12 | 220 | 8.55 | 414.46 |

All three input files are byte-identical to the corresponding reference;
the runner explicitly selects CIAH and leaves the deck unchanged. The 241
case is the qualification guess (`5dbd9c72...`), not `241_vb.xmi`'s different
guess. Its final energy differs by `9.2e-11 Eh`, HVP falls 34.1%, AH Ritz
solves fall 78 to 46, and evaluated keyframes stay at 13. Process wall falls
2.03 to 1.47 s in these single Hanhai25 observations. Both nodes use EPYC
9755, but node/load and compiler differ (reference GCC 14.1 versus 14.3).

The 240/FeCl2 rows are older-version context, not an isolated ablation of
this patch: intervening commits also changed the transported finite model
and contraction kernels. Hanhai22 timings and different core counts must
not be used to attribute a speedup to this proposal. The current 240 RSS is
also higher than that old reference; no claim of a memory reduction follows
from these records. The path itself adds no owned direction arrays.

The first implementation owned two additional direction arrays and was
tested in jobs 291706 and 291711; logs remain in `native/`. Borrowing their
already existing images preserves all five step/HVP counts and final
energies within `2e-12 Eh`. Time variation between those runs is not credited
as a separate optimization.

## Remaining target

FeCl2's last three macros use 25, 27 and 29 HVPs (81 of 170). MnF2's last
four use 18, 20, 22 and 26 (86 of 160). Their actual microstep logs show no
clipping or replay in those tails. Thus segmented advancement cannot remove
that remaining cost: a subsequent improvement must address how the AH model
and its complement are used in the unclipped regime, without changing the
outer accuracy or simply tightening the inner Ritz solve. No iron-porphine
long run was launched in this round.

## Reproduction

Use the existing `../ciah-subspace-20261017/build.sbatch` on the qualified
execution mirror. Place the frozen decks under this directory's remote
`inputs/{240,241,FeCl2,F2,MnF2}.xmi`, then submit `run.sbatch SOURCE_DIR` with
`--array=0-4` from this directory. The first three decks are preserved in
`results/diagnostics/vb-adapters-20261016/{240,241,FeCl2}/input.xmi`;
F2 and MnF2 are repository testdata. The CSV records their full hashes and
artifact locations. Do not rerun the reference binaries.
