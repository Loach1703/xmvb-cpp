# Oriented factor and residual-generator adjoints

## Implemented kernel and production boundary

`mcscf/ci/factor_generator.cpp` implements two exact contractions:

- `factor_generator_adjoint`: `<X, G(kappa)Y>`;
- `factor_residual_generator_adjoint`: `<(H-E)X, G(kappa)Y>`.

Only weighted single-spin panels `X F_right`, `X^T L_right` and their
residual-weighted counterparts are formed. Products with X use source-pair
moments. The Hamiltonian case reuses the existing `FactorActionWorkspace`
on weighted and ket panels, then reuses its consumed raw-channel buffers for
the crossed moment/gather. No new sigma implementation, dense product array,
cofactor path or per-thread spin/product array is added.

`DirectCiSigmaAction::spin_generator_adjoint` retains oriented creation and
annihilation indices. Ordinary left/right interchange requires a transpose;
Hamiltonian-left interchange is not permitted, since H and an oriented
generator generally do not commute. Tests explicitly use different left/right
ranks, nonsymmetric generators and unequal spin spaces.

Accepted-state transition densities now use this factor kernel when its own
moment/link work estimate is smaller. They do not reuse the unrelated
Hamiltonian-channel cost model. Selected-state coefficient responses and
off-diagonal state-average multipliers keep their original meanings. No
optimizer strategy, coordinate, physical space or tolerance changed.

The residual kernel is an implemented and validated prerequisite for removing
accepted full residuals. It is NOT yet substituted into those consumers that
already have a cached residual: doing so could repeat Hamiltonian work.
Full accepted states/residuals and their remaining HVP consumers still need
migration. This round does not establish bounded memory for the whole SCF.

Theory, proof, leading work and numerical-workspace bounds:
Section8.3 of `article/source_column_hamiltonian.tex` and its compiled PDF.

## Verification

Hanhai22 anode26,8 bound cores for kernels, GCC14.1.0/CMake3.30.5, Release
`-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`, BLAS1. Final numerical/SCF
job1260854 passed9/9 CTests. Additional strided-view pairing checks passed
9/9 in1260871;1260873 also passed9/9 before its larger native test.

Coverage includes unequal/empty beta sectors, different bra/ket ranks,
full and very small source/channel workspaces, directional providers,
non-owning row-strided views and independent generator/adjoint pairing.
The multistate test compares every retained off-diagonal transition density
against a complete-product oracle for equal and unequal weights; existing
joint orbital/structure HVP differences and gauge-multiplier checks still pass.

All builds and numerical execution used Slurm compute allocations. No old
executable was rerun. Measurements use native clocks and GNU time, not
unavailable profiler sampling or inferred hotspot percentages.

## Kernel measurements

One physical covalent-structure algebraic fixture per active space. These are
not molecular SCF or complete orbital-HVP timings. The following factor and
reference times share the same reference process. Factor times include source
generation/exterior transforms. Reference amplitude preparation is excluded.

| CAS | Bra/ket factors | Transition reference/factor s | Gain | Residual generation+contraction reference/factor s | Gain |
|---|---:|---:|---:|---:|---:|
|12/12|64/32|0.007318/0.002310|3.17x|0.112053/0.064527|1.74x|
|14/14|128/64|0.145981/0.018431|7.92x|2.018170/0.619571|3.26x|

Generator error norms are4.52e-14/3.22e-13; residual-generator error norms
are1.59e-14/2.00e-13. The returned orbital matrices are not artificially
symmetrized. Native-only process RSS is53.78/61.66 MiB. The CAS14 residual
workspace bound is60,739,872 bytes, below the unchanged64 MiB budget,
including provider scratch and the small thread-local orbital reductions.

The important cached-residual distinction is measured separately:

| CAS | Sigma formation s | Contraction of already available residual s | Factor residual total s |
|---|---:|---:|---:|
|12/12|0.090703|0.021350|0.064527|
|14/14|1.288312|0.729858|0.619571|

CAS12's cached-residual contraction is cheaper than recomputing the factor
Hamiltonian action. This is why the new residual kernel must accompany removal
of the old residual producer/consumers, not be bolted onto that cached path.
No claim of a1.74x/3.26x complete-SCF speedup is made.

## Complete 240 regression

Job1260854 also ran COV,ION(0:1),ION(0:2),ION(0:3),FULL to convergence.
Same node/eight cores and exact-libcint NEO/Davidson settings as saved1260422:
projected gradient1e-3, energy change1e-7 Eh, Hessian-diagonal preconditioning.

| Space | Steps | HVP | Saved/current SCF s | Saved/current RSS MiB | Final energy Eh |
|---|---:|---:|---:|---:|---:|
|COV|3|20|6.676/6.491|1528.6/1518.4|-343.239211501260|
|ION(0:1)|3|27|8.761/8.821|1517.8/1508.0|-343.407736636064|
|ION(0:2)|4|32|10.834/10.673|1518.7/1546.5|-343.443415326852|
|ION(0:3)|4|32|10.983/11.059|1515.6/1522.4|-343.446256564846|
|FULL|4|30|9.623/9.465|1525.7/1522.4|-343.446302545201|

Iterations/HVP counts are unchanged; maximum final-energy difference is
1.99e-12 Eh. Whole-task time/RSS remain essentially unchanged for these
broad-support cases in single native runs. Old results come from CSV.

## Larger native qualification and next work

Job1260873 submits a16/16 single-covalent native kernel, with128 representative
ket sources,256 physical bra sources and12870 strings per spin. It does not
construct a full12870-by-12870 reference. Its timing/RSS must be collected
separately;12/14 timings are not evidence of16/16 performance.
At the64 MiB budget the present planner selects source widths1/1 and channel
width105. Its two channel panels imply131072 bra-source loads. This is a
code/operation-count diagnosis, not a measured hotspot percentage. A candidate
with full256-source bra residency and channel width20 fits the same budget;
joint planning must weigh the avoided source generation/output updates against
the extra channel traversals. The running native job is not restarted.

Next required work is the complete accepted-state/HVP consumer migration,
including fused integral/generator moments and exact structure images without
full states/residuals. The64 MiB source/channel planner must also be qualified
at the larger scale: repeated source regeneration can dominate when a large
channel panel leaves little room for bra sources. Do not fix this by silently
enlarging the budget or treating the selected subspace as closed.

Raw evidence and CSV: `factor_generator.csv`, `factor_generator_240.csv`,
`factor-generator-1260854`, `factor-generator-1260871`, associated application
logs, hashes and Slurm accounting. The parser is
`tests/mcscf/summarize_factor_gradient.py`.
