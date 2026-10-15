# Directional spin factors in the integral part of an HVP

## Mathematical and code change

The orbital representation derivative is the exact factor sum
`dY = (G_alpha(kappa)L) F^T + L (G_beta(kappa)F)^T + T_alpha W dC T_beta^T`,
where `kappa = dR R^-1`. This differentiates the CI representation, not the
physical finite nonorthogonal orbital update. Original selected VB support
need not be closed under either generator. No intermediate projection,
amplitude cutoff, approximate rank or convergence threshold is introduced.

`mcscf/ci/generator_factors` loads these orbital terms on bounded source panels;
`sum_factors` adds the selected-structure coefficient response. The common
single-spin generator preserves oriented creation/annihilation indices, since
kappa is not generally symmetric. The common integral-adjoint kernel now
accepts these sources directly. Its cost comparison includes all directional
source families and the source-pair moments, rather than assuming the accepted
and directional ranks are equal.

`StructureAction::active_integral_adjoint_direction()` uses this order when its
arithmetic estimate is smaller. The common support expansion also replaces
the duplicate original-support expansion loop in the forward source consumer.
Loader scratch is part of the same numerical workspace budget, not hidden
outside it. Exact product and factor orders still evaluate the same operator;
the planner does not admit/reject a physical selected space based on memory.

Theory and workspace bounds are in Section8.2 of
`article/source_column_hamiltonian.tex` and its compiled PDF.

## Avoiding repeated directional source generation

The first CAS14 directional-factor order was only slightly faster than the
complete-product reference. Balanced source panel selection left the bra
nonresident, regenerating its exterior maps and generator images for each ket
panel. This is unnecessary arithmetic, not insufficient parallelism.

If the complete bra and one ket column fit the unchanged workspace budget,
reserve the bra first and maximize the ket width in the remainder. The leading
Gram/moment arithmetic is independent of this width; bra generation then occurs
once. If it cannot fit, the same equations stream both source families.
For CAS14 the resulting panel dimensions are bra256/ket3/channels105, with
56,655,360 bytes of explicit slab/loading workspace, below64 MiB. Phase clocks
show bra/source loading about0.035 s, scatter0.118 s, covectors0.356 s and
Gram0.119 s. These are explicit stage clocks, not sampled profiler percentages.

## Validation

Final job **1260805**, Hanhai22 **anode26**, eight bound cores,
GCC14.1.0/CMake3.30.5, production Release `-O3 -DNDEBUG -g1
-fno-omit-frame-pointer`, BLAS1. All **9/9 CTests passed**.

New checks include nonsymmetric generators, unequal alpha/beta spaces, empty
beta occupation, panels crossing source-family boundaries, joint source/channel
streaming and loading-scratch accounting. Production HVP checks compare the
complete gradient derivative against central differences for orbital-only and
combined orbital/structure directions. A sparse unequal-spin selected space
exercises the new order explicitly. A two-state selected-space test compares
the full HVP pullback with an independent complete-product contraction for
equal and unequal state weights and nonzero off-diagonal gauge multipliers.
Existing parity, mixed-support, metric, coupling and transport tests also pass.

## Native directional integral-adjoint kernels

These are one-covalent-structure algebraic CAS fixtures, not total molecular
HVP or SCF times. Native RSS uses separate reference-free processes. The
reference column evaluates the complete-product integral adjoint only; amplitude
and directional-carrier preparation are excluded from its timer. The factor
timer includes source/exterior/generator loading.

| CAS | Strings/spin | Direction/ket factors | Native warm s | Reference-process factor s | Complete-product reference s | Native RSS MiB | Integral error norm |
|---|---:|---:|---:|---:|---:|---:|---:|
|12/12|924|128/32|0.0666|0.0610|0.1287|74.91|7.90e-17|
|14/14|3432|256/64|0.6364|0.6279|1.1812|75.22|1.09e-16|

Within each reference process, the factor order is **2.11x / 1.88x** faster.
The preceding nonresident candidate1260782 took0.0873/1.1026 s in native-only
processes, with74.34/83.35 MiB RSS. The new native order is1.31x/1.73x faster
at the same budget; CAS14 RSS falls by about9.8%. Raw observations from the
intermediate1260788 job are also retained in `factor_direction.csv`, rather
than selecting its best timing for the final table. No historical executable
was rerun. No claim of universal wall-time speedup follows from these cases.

## Complete 240 convergence regression

The final1260805 allocation also builds the application and runs the same five
240 subspaces to convergence. Settings match saved1260422: anode26/eight cores,
exact libcint, NEO/Davidson, Hessian-diagonal preconditioner, projected-gradient
threshold1e-3 and energy-change threshold1e-7 Eh. Old data are read from CSV.

| Space | Steps | HVP | Saved/current SCF s | Saved/current RSS MiB | Final energy Eh |
|---|---:|---:|---:|---:|---:|
|COV|3|20|6.676/6.477|1528.6/1534.8|-343.239211501261|
|ION(0:1)|3|27|8.761/8.634|1517.8/1518.9|-343.407736636064|
|ION(0:2)|4|32|10.834/10.840|1518.7/1521.3|-343.443415326852|
|ION(0:3)|4|32|10.983/10.976|1515.6/1498.5|-343.446256564845|
|FULL|4|30|9.623/9.502|1525.7/1523.1|-343.446302545203|

Iterations/HVP counts are unchanged; maximum final-energy difference is2.05e-12
Eh. These broad-support cases have essentially unchanged total time/RSS in
single runs. They verify nonregression, not a complete-SCF factor speedup.

## Still incomplete

The new source order removes the full-product Gram from the HVP's integral
adjoint when advantageous. Full orthogonal states, sigma/residual arrays and
their generator-gradient derivatives remain in the surrounding HVP pipeline.
The selected-state response producer also remains dense. Their remaining
allocations cannot be excused by this kernel's bounded workspace.

Next milestone remains the unsymmetrized factor/generator VJP and residual
contraction, then removal of these obsolete full carriers after all consumers
are migrated. See repository-root `TODO_ORTHOGONAL_CI.md`.

Data: `factor_direction.csv`, `factor_direction_240.csv`, raw logs and hashes.
Parser: `tests/mcscf/summarize_factor_gradient.py`. Timings and RSS are native;
the earlier unavailable gprofng capture is not represented as profiler evidence.
