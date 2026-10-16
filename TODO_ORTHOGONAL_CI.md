# Exact orthogonal-CI contraction for arbitrary VB structures

Goal: preserve the selected VB variational space and structure-space Davidson,
while evaluating the complete forward/backward chain by Slater--Condon spin
operations and bounded workspaces. No selected-space closure, amplitude cutoff,
or memory-based admission to a different physical representation is allowed.

## Current audit

- All three active-space forward entries now complete the computational
  one-spin topology before constructing `StructureAction`. Previously only
  synthetic selected-space tests called that constructor; real incomplete
  string lists could therefore miss the new path. Original input expansion
  maps and the selected VB variational space remain unchanged. The completion
  stores one-spin metadata, not all spin pairs or product amplitudes.
- Fixed-point H/M actions have two exact contraction orders: generated product
  tiles and two-sided source factors. Neither needs a dense transformed product
  array. The arithmetic planner chooses their order, not the physical space.
- `prepare_active_adjoint()` now retains raw coefficients in structure space,
  not as a product matrix. Its integral adjoint can stream spin factors using
  the exact factor/product arithmetic planner. Source-order accepted states
  now omit orthogonal/residual product arrays and use factor derivatives.
  Broad product-order states still retain those arrays.
  Transition one-body densities now have an exact streamed factor order,
  chosen by their own moment/link work rather than the Hamiltonian-channel cost.
- Product-order `apply_integral_direction()` allocates stacked product arrays for the base
  coefficients, sigma, coefficient direction and overlap direction.
- Product-order `active_integral_adjoint_direction()` still consumes
  those dense carriers for the residual/generator terms. It now uses
  exact directional source factors for its integral adjoint when their work is
  smaller. Tiling only the Davidson H/M action is therefore not a
  bounded-memory VBSCF implementation.
- `SelectedStateDeterminantCoefficients::coefficient_matrix` is itself dense;
  changing only its downstream consumers cannot remove this allocation.
  The accepted-point producer now skips those arrays for direct CI; the NEO
  response-work estimate uses topology/state dimensions, independent of numerical
  array residency. `test_keyframe --accepted-state16` exercises the complete
  gradient producer and confirms zero raw/orthogonal/residual product doubles
  at 12870 strings per spin (117.48 MiB process peak, constructed 16-AO input).
  The direct-CI response interface now also borrows structure coefficient columns:
  `apply_structure_response_adjoint()` and the gauge/relaxed-response branches
  of `exact/apply.cpp` no longer create a raw selected-state product bundle.
  Source-order consumers expand only original touched support. Product-order
  consumers generate response columns and consume their sigma/generator adjoints
  in the same tile sweep, with no full response carrier. Accepted arrays remain
  resident on this broad-support path. Symbolic-support column panels share row
  indices across RHSs; their bounded caches retain zeros inside the support but
  do not store the forbidden rows.
  Pure response adjoints omit unused structure H/M output projection.
- Without a supplied Davidson diagonal, `StructureAction` still evaluates a
  double loop over the expansion terms of each structure through spin-pair
  matrix elements. This must not survive as an accidental cofactor bottleneck
  after H/M itself has moved to orthogonal CI. Trace precomputed-diagonal
  producers too; do not infer their runtime share from this code inspection.
- `RestrictedMetricInverse` also enumerates the full product carrier for
  occupation/complement metadata and reserves full-carrier RHS workspace.
  Audit and migrate this consumer as well; H/M tiling alone does not bound it.

## Ordered milestones

Production correctness gate: job1261206 captured an actual LOFLEA CAS(12,12)
metric-symmetry failure and its exact raw orbitals after13 accepted steps.
The metric norm used a physical tangent map while its action used independent
projector simplifications. The action now uses the strict adjoint of that same
map, with no changed certificate, coordinate whitening or new inverse solve.
At the identical point, the action/feature-Gram gap drops from1.14e-9 to2.56e-11
and skew from1.90e-10 to8.55e-12 (jobs1261210/1261230). New cancelling-representative
and accepted-cache tests pass. A bounded original-input capture reaches25 steps
without that error but is not converged (job1261218); do not report complete
CAS12 molecular qualification. Keep its slow convergence distinct from the
qualified source/product derivative kernels and the captured-point repair.
See `article/orbital_metric_pullback.tex` and the metric consistency report.

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
   Generated product tiles now also have one shared sigma/integral-VJP graph:
   scatter each input once, form the channel Gram, reuse its right slab for K*D,
   and gather directly into a consumer. A second generated input reuses the
   consumed input halo, not another product allocation. Both exchange parities
   and unequal spin spaces are covered. The dense-carrier wrapping experiment
   did not improve production convergence/performance and is removed; migrating
   generated readers and residual/generator consumers remains necessary.
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
   panels without a residual carrier, reusing the shared channel sweep. Source
   states now remove accepted orthogonal/residual storage; broad product-order
   states still need tile-consumer migration before deleting their storage.
   Route exact structure diagonals through an appropriate contraction consumer,
   without probing every global basis vector or adding quadratic structure work.
5. [ ] Differentiate that same factor graph for orbital HVP and CI responses.
   Include factor tangents, integral tangents, exterior-transform pullbacks and
   their directional terms. Arbitrary selected spaces may not be treated as
   invariant. Do not add relaxed CI response to a coupled action a second time.
   Partial: `generator_factors` represents both oriented orbital-generator
   terms. `SpinSources::tangent_source` merges coefficient changes into their
   non-anchor spin factors, so combined directions use at most2r rather than3r.
   HVP integral adjoints consume this exact source sum; its temporary loading
   panel is included in the numerical workspace budget. Source-state orbital
   H/M directions and complete gradient/HVP pullbacks are now factor-based,
   including left/right representation changes, delta-H and gauge terms.
   The H-image sweep also collects the small generator-gradient direction;
   subsequent HVP contraction reuses it instead of repeating Hamiltonian work.
   Broad product-order directional consumers still need migration.
   The coupled structure-coefficient response now has a bounded factor order:
   expand its original support once, reuse the common H/M source sweep, and
   collect both generator terms from the shared Hq/HY images without a residual
   Hamiltonian replay; only Hq is projected to structures. This
   consumer needs no accepted or response product array; broad-support cases
   retain the exact product order when its total adjoint/image work is lower.
   The exact-HVP producer now passes both horizontal and normalization/gauge
   coefficient responses directly in structure space. The old raw-response
   public overload is removed, and the source-only adjoint uses the same sweep
   with structure output collection disabled. Nine focused tests passed in
   job1261329; complete 240 subspaces and arbitrary 1/3-structure SCF regressions
   passed in job1261339 without changing steps or HVP counts.
   Production CAS16/16 probe1261340 now covers the accepted-gradient producer,
   one complete orbital/coupling HVP (including normalization response), and
   central differences of two relaxed gradients in the same raw linear lift.
   It passes at1.73e-9 absolute HVP error, with155.70 MiB process peak and no
   accepted product arrays or raw response rebuild. This is a constructed
   16-AO Hamiltonian, not a converged molecular SCF benchmark. Its HVP still
   costs153.49 s; stage attribution and contraction reduction remain necessary.
   Jobs1261360/1261367 passed mixed row/column anchor reconstruction, unequal
   spins, both exchange parities, generator-only/coefficient-only/combined
   tangents, cross-panel boundaries, factor/product derivatives and9 complete
   SCF regressions. Production probe1261368 retains1.73e-9 HVP difference error
   and155.94 MiB RSS, but costs154.54 s versus saved153.49 s: the exact2r
   simplification is not an end-to-end speedup. Its55.01 s selected direction
   and99.53 s backward phase identify where further work is needed.
   The equivalent parity-placement order is now implemented. Existing phase
   capture1261387 identified106496 planned directional-bra column loads and
   52.96 s provider generation. New job1261411 keeps the accepted physical bra
   resident under the same64 MiB cap:256 planned loads and0.65 s generation.
   The complete production HVP drops from155.80 to100.29 s, with155.39 MiB RSS
   and unchanged1.73e-9 difference error. No oriented residual adjoint is swapped.
   Jobs1261395/1261413 pass9 focused tests and9 full240 regressions.
   FeCl2 control1261419 reaches the correct energy, but needs16/271 steps/HVP
   versus12/211 at older metric checkpoint1261264. Do not claim universal
   end-to-end improvement; this unequal-spin control has no common exchange
   quotient, so keep its trajectory change separate from the new parity order.
   Next derive direct accepted-adjoint reuse for the known within-state
   normalization response, using its upstream algebraic frame rather than a
   numerical amplitude threshold. Broad product-order streaming remains open.
   For a single root with known Q_g=aC, the coefficient-only integral and
   residual-generator adjoint derivatives equal2a times their accepted values.
   Reuse the accepted active-integral adjoints, not the total orbital gradient
   (which also contains inactive-only terms). Keep the gauge multiplier and
   orbital-map derivatives separate. Obtain a from the core coupling's small
   selected-frame solve, never by testing whether numerical vectors are parallel.
   Weighted multistate totals alone do not contain all cross-state adjoints.
   Completed: retain the upstream small gauge frame instead of eagerly forming
   its full coefficient response. Single-root direct-CI HVPs reuse exactly 2a
   times the accepted active-integral adjoints. Multiplier/map derivatives and
   general multistate frames remain. Jobs1261467/1261483 pass11 focused tests,
   9 complete 240 optimizations and scalar/batch fixed-CI HVP differences with
   inactive orbitals and1/2 roots. Native CAS16/16 probe1261462 reduces HVP
   100.29 ->76.26 s versus saved1261411, at155.11 MiB peak RSS and1.73e-9
   difference error, removing both extra normalization residual-generator
   sweeps. FeCl2 control1261472 reaches-2181.624464864478 Eh in13 steps/235 HVP;
   the older12-step/211-HVP correct checkpoint remains a performance target.
   See `normalization_reuse.md` and the automatically parsed CSV tables.
6. [ ] Remove obsolete full-product production consumers after migration.
   Explicit full-carrier arrays remain only as independent small-system test
   oracles; broad-support production uses generated product tiles, not a dense
   full-product allocation.
   Next priority: give the existing product-tile sigma/adjoint graph generated
   coefficient readers, so broad-support accepted gradients and orbital
   directions stream too. Share the H-image consumer, retain oriented generator
   terms and their support halo, and do not instantiate another full carrier.
   Structure-response producer migration is complete for direct CI. The separate
   non-direct pair consumers retain their own explicit numerical representation.
   Cholesky directions preserve the generated downset S: delta-R R^-1 is
   upper triangular, so both orbital and selected-coefficient dY remain in S.
   The physical overlap pullback consumes only upper generator entries; those
   entries and their derivatives depend only on residuals projected to S.
   Verify and exploit this certificate in the streaming consumer, not a blanket
   projection of general oriented generators or of arbitrary selected VB spaces.
   Completed for coefficient response: H/M, sigma/adjoint and physical generator
   pullback share generated panels and reuse accepted Y/R. Compact exact-support
   rows remove padded-zero cache thrashing without increasing cache budgets.
   Hanhai22 job1261729 validates a CAS14/14 selected132 response against an
   independent complete-carrier differential (image1.11e-15, adjoint6.45e-16).
   Cold response falls293.93 ->3.33 s versus saved1261623, process RSS1661.66
   ->438.14 MiB; all11 focused tests pass. This is a coefficient-response kernel
   result, not a complete-SCF speedup. See `product_response.md` and its CSV.
   Accepted broad-support Y/R and orbital-direction carriers remain the next
   ownership migration. The single-source factor path's older70 ms record also
   remains a separate target; this change does not improve every contraction order.
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

Resolved large-scale scheduling defect: saved job1260873 used CAS16/16 with256
physical bra and128 ket factors and selected widths1/1/105, causing131072 bra
source-column loads. Its residual kernel cost174.37052 s with118.80 MiB RSS.
The joint scalar-work/traffic planner in `mcscf/ci/factor_plan.cpp` now selects
widths256/1/20 under the same64 MiB budget. Job1260948 measured12.12480 s,
118.52 MiB RSS,256 loaded columns in one provider call. The channel-sorted
link traversal uses no additional offset table or per-thread large arrays.
Its16/16 low-memory directional integral/number identity passed at8.33e-17;
no complete product matrix was formed, so full-matrix error remains unmeasured.
Use these saved results rather than rerunning historical candidates. Remaining
work is the accepted-state/HVP consumer migration and fused integral/generator
moments; the residual consumer still restricts partial-channel ket width to one.

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
