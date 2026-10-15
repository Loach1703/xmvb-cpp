# Exact orthogonal-CI contraction for arbitrary VB structures

Goal: preserve the selected VB variational space and structure-space Davidson,
while evaluating the complete forward/backward chain by Slater--Condon spin
operations and bounded workspaces. No selected-space closure, amplitude cutoff,
or memory-based admission to a different physical representation is allowed.

## Current audit

- Fixed-point H/M actions have two exact contraction orders: generated product
  tiles and two-sided source factors. Neither needs a dense transformed product
  array. The arithmetic planner chooses their order, not the physical space.
- `prepare_active_adjoint()` now retains raw coefficients in structure space,
  not as a product matrix. Its integral adjoint can stream spin factors using
  the exact factor/product arithmetic planner. Orthogonal coefficient matrices
  and residuals remain on the full product carrier for every selected state.
  Transition one-body densities now have an exact streamed factor order,
  chosen by their own moment/link work rather than the Hamiltonian-channel cost.
- `apply_integral_direction()` allocates stacked product arrays for the base
  coefficients, sigma, coefficient direction and overlap direction.
- `active_integral_adjoint_direction()` and `contract_response()` still consume
  those dense carriers for the residual/generator terms. The former now uses
  exact directional source factors for its integral adjoint when their work is
  smaller. Tiling only the Davidson H/M action is therefore not a
  bounded-memory VBSCF implementation.
- `SelectedStateDeterminantCoefficients::coefficient_matrix` is itself dense;
  changing only its downstream consumers cannot remove this allocation.
- Without a supplied Davidson diagonal, `StructureAction` still evaluates a
  double loop over the expansion terms of each structure through spin-pair
  matrix elements. This must not survive as an accidental cofactor bottleneck
  after H/M itself has moved to orthogonal CI. Trace precomputed-diagonal
  producers too; do not infer their runtime share from this code inspection.
- `RestrictedMetricInverse` also enumerates the full product carrier for
  occupation/complement metadata and reserves full-carrier RHS workspace.
  Audit and migrate this consumer as well; H/M tiling alone does not bound it.

## Ordered milestones

Target ownership, without a new generic SCF class hierarchy:

| Layer | Owns | Must not own |
|---|---|---|
| VB support plan | structure incidences, touched products, minimum source cover, causal transforms | a transformed product array |
| Accepted numerical state | structure/state coefficients, orbital/integral data, bounded source residency | one dense product copy per state/response |
| Shared `mcscf/ci` factor sweep | spin scatter, channel mixing, gather, reusable tile workspace | VB localization or optimizer policy |
| Consumers | H/M sampling, integral/factor adjoints, differentiated contractions | separately reimplemented spin kernels |
| Optimizer | structure-space Davidson and existing orbital coordinates | carrier representation decisions or duplicate CI response |

The shared numerical module is `mcscf/ci/factors`, not another parallel
implementation of the source loop. The forward sweep and integral adjoint are
implemented there; factor adjoints and tangent consumers remain next. A source
provider generates panels from sparse coefficients; it must not first build
full `U * r` factors when a broad support makes those product-sized.

1. [x] Replace beta-only sources by an exact minimum row/column support cover.
   Assign each original input edge to one factor; preserve the complete original
   output support. A maximum matching certifies the minimum factor count among
   support-only coordinate-anchored decompositions. No numerical rank estimate.
2. [x] Keep causal exterior-adjoint plans with the immutable accepted-point
   support. Do not rebuild their dependency graphs for every Davidson vector.
   Structure incidences now store the support-group index directly, removing
   repeated binary searches in source/product/overlap back-projection without
   adding an index array or enlarging each incidence.
3. [x] Extract one shared spin-factor channel sweep from H/M into `mcscf/ci`.
   Reuse bounded backing allocations across right-hand sides. Add exact streamed
   integral adjoints, transferring self-adjoint generators onto the smaller ket
   source family to avoid repeated bra scatter. Both consumers reuse the common
   spin kernels and Hamiltonian-kernel pullback; factor VJP remains in item 4.
4. [ ] Introduce sparse selected-state support coefficients at their producer.
   Replace dense accepted raw/orthogonal/residual carriers by a shared support
   plan, numerical source factors and streamed images. Migrate accepted integral
   gradients and transition one-body adjoints together, including unequal spins
   and state-average cross-state gauge multipliers.
   Partial: accepted raw state caches now hold structure coefficients. The
   direct-CI gradient entry no longer builds an unused raw selected-state bundle.
   Integral gradients choose the exact factor or product contraction by work;
   their nonstationary full-carrier generator term is retained, not discarded.
   `factor_generator_adjoint` now supplies oriented transition densities to this
   producer. `factor_residual_generator_adjoint` computes residual-weighted
   panels without a residual carrier, reusing the shared channel sweep; remaining
   residual consumers must migrate before deleting accepted residual storage.
   Route exact structure diagonals through an appropriate contraction consumer,
   without probing every global basis vector or adding quadratic structure work.
5. [ ] Differentiate that same factor graph for orbital HVP and CI responses.
   Include factor tangents, integral tangents, exterior-transform pullbacks and
   their directional terms. Arbitrary selected spaces may not be treated as
   invariant. Do not add relaxed CI response to a coupled action a second time.
   Partial: `generator_factors` represents both oriented orbital-generator
   terms, and `sum_factors` includes selected-structure coefficient responses.
   HVP integral adjoints consume this exact source sum; its temporary loading
   panel is included in the numerical workspace budget. Full carrier residual
   and generator-gradient derivatives still need migration.
6. [ ] Remove obsolete full-product production consumers after migration.
   Explicit full-carrier arrays remain only as independent small-system test
   oracles; broad-support production uses generated product tiles, not a dense
   full-product allocation.
7. [ ] Validate the complete SCF chain on 240 COV/ION/FULL and selected CAS(12,12)
   spaces, then larger covalent cases. Report convergence, energy, HVP count,
   total wall time and native peak RSS, not only a fixed-point kernel speedup.

The source cover is generated once from touched input edges, never by scanning
all alpha/beta products. Source factors and channel panels are numerical data;
the selected structure incidence, minimum cover, spin links, and causal
transform plans are immutable shared metadata. Numerical workspaces belong to
one evaluation, not to individual OpenMP threads.

For derivative migration, bilinear left/right factor contractions generally
introduce source-pair Gram matrices of size `r_left * r_right`. Their work is
not automatically linear in the number of factors. An exact cost planner must
compare this order with generated product tiles for each consumer. Both are
evaluation orders of the same physical operator, not fallback approximations.
Factor tangents obey `dY = dL F^T + L dF^T`; finite orbital updates and the
selected-space pullback must remain consistent with that differentiated graph.
Hamiltonian channels are symmetric packed one-body generators, but orbital
representation pullbacks also need the unsymmetrized generator adjoint. Retain
the oriented creation/annihilation indices from the existing spin links; do not
silently symmetrize that adjoint when sharing the channel implementation.

## Required numerical and resource evidence

Immediate larger-scale qualification: native generator job1260873 is CAS16/16
with256 physical bra and128 ket factors. The present residual planner chooses
`left_tile=1,right_tile=1,channel_tile=105` at64 MiB. Two channel panels cause
`128 * (2+2) * 256 = 131072` bra-source loads. A feasible resident-bra point
exists at `left_tile=256,right_tile=1,channel_tile=20`, with the same budget.
Do not assume the larger channel tile is optimal. Joint source/channel planning
must count exterior/generator regeneration, repeated spin-link traversals and
repeated output updates across bra panels. Preserve the running job and collect
its actual time/RSS; do not rerun historical binaries or increase the budget.

- Exact source reconstruction, both spin-exchange parities, mixed parity,
  unequal spin lists, ionic diagonal support and multiple right-hand sides.
- H/M images versus an independent complete-carrier oracle; integral-adjoint
  directional identities; orbital gradient differences; HVP differences and
  bilinear symmetry; equal-weight multistate gauge tests.
- Run builds/tests only in Hanhai22 Slurm allocations. Reuse saved historical
  results rather than rerunning old binaries. Record input/source/binary hashes,
  affinity, native time/RSS and full allocation ownership.
- A source-factor improvement is not automatically a gradient/HVP improvement.
  Retain this distinction until milestones 4 and 5 are complete.

Theory: [TeX derivation](article/source_column_hamiltonian.tex),
[PDF](article/source_column_hamiltonian.pdf).
