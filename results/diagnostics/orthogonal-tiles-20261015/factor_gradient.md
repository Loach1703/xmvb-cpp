# Shared factor contractions and accepted-state gradient integration

## Implemented boundary

The fixed-point H/M factor sweep is now in `mcscf/ci/factors`, shared independently
of VB data types. VB supplies bounded source panels and sparse sampling/projection
callbacks. The numerical workspace is allocated once and reused across RHS;
spin/channel scatter, GEMM and gather are no longer separately implemented in
the VB source consumer.

The new exact integral adjoint contracts two spin-factor families without a
complete product amplitude. Self-adjoint generators are transferred to the
smaller ket family. Weighted bra bases are accumulated once per ket panel;
channel covectors consume all bra panels before the Q-channel Gram. This
eliminates repeated bra scatters and repeated leading channel contractions.
The explicit numerical slab budget stays 64 MiB; no tolerance, localization,
structure-space or spin restriction changes.

Accepted-state preparation chooses the exact factor or product integral
contraction by an arithmetic estimate including the source-pair term. It now
retains raw coefficients in structure space, not as another dense carrier.
The direct-CI gradient entry no longer produces an unused raw selected-state
bundle. Full carrier residual terms in the exterior-map gradient are retained:
stationarity in an arbitrary selected space is not full-carrier stationarity.

Standalone product sigma and fused sigma/adjoint now use one shared channel
sweep. Certified common exchange parity is propagated to derivative sigma
calls. A parity call returns the full image, not a representative image. The
duplicated sigma implementation and its thread-private beta reductions were
deleted, not kept as a fallback. No new copy per thread is introduced.

## Validation and environment

Final unit/kernel job **1260725**, Hanhai22 **anode26**, eight bound cores,
GCC14.1.0/CMake3.30.5, Release `-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`,
BLAS1. All **9/9 CTests passed**. Preceding integration job1260678 also passed9/9;
1260725 additionally checks the accepted-state directional gradient/HVP after
raw-cache removal and sigma unification.

Tests cover unequal spins, mixed source covers, both exchange parities, partial
source/channel panels, zero source rank, different left/right factors, non-owning
views, workspace reuse, and an independent integral-direction identity. A
selected CAS(8,8) gradient fixture has a nonzero full-carrier residual and checks
the complete integral/exterior pullback plus its central-difference HVP.

These are native stage clocks and GNU time records, not profiler percentages.
The saved compute-node gprofng probe failed; no sampling attribution is claimed.

## Physical covalent-state integral-adjoint kernel

One covalent structure, complete one-spin carrier, nonorthogonal active orbital
fixture; not a molecular SCF calculation. Final1260725 measurements below use
separate native-only processes for RSS. The full-product reference is the
current independent contraction order, not a rerun of a historical executable.
Reference time excludes amplitude preparation; factor time includes source
generation and transformation.

| CAS | Strings/spin | Bra/ket factors | Native warm s | Full-product reference s | Same-process factor s | Native RSS MiB | Integral error norm |
|---|---:|---:|---:|---:|---:|---:|---:|
|12/12|924|64/32|0.0565|0.1016|0.0664|74.64|3.90e-15|
|14/14|3432|128/64|0.4413|1.1199|0.4384|82.52|8.99e-15|

Within the reference processes the factor order is **1.53x / 2.55x** faster.
The CAS14 explicit factor workspace is65,773,056 bytes, below64 MiB. Process
RSS also includes spin/exterior metadata, integral kernels and returned data.

The implementation alternatives are retained in `factor_adjoint.csv`, not in
production code. CAS14's initial bra-scatter order cost4.657 s, a first batched
layout cost3.359 s; the final generator-transfer order costs0.441 s. Thus merely
batching the first formula was insufficient. The unused alternate scatter
layout was removed. No old executable was rerun to obtain this comparison.

## Complete 240 convergence regression

Job **1260738**, same anode26/eight-core configuration as saved baseline1260422,
same exact libcint input, NEO/Davidson, Hessian-diagonal preconditioner,
projected gradient1e-3 and energy change1e-7 Eh. Each current subspace ran once
to convergence. Old results are read from `source_240.csv`, not recomputed.
Times below are SCF-loop wall times; RSS is whole native process peak.

| Space | Structures | Steps, old/new | HVP, old/new | Saved/current SCF s | Saved/current RSS MiB | Final energy Eh |
|---|---:|---:|---:|---:|---:|---:|
|COV|14|3/3|20/20|6.676/6.545|1528.6/1518.7|-343.239211501260|
|ION(0:1)|294|3/3|27/27|8.761/8.799|1517.8/1527.3|-343.407736636063|
|ION(0:2)|1134|4/4|32/32|10.834/10.702|1518.7/1519.3|-343.443415326850|
|ION(0:3)|1694|4/4|32/32|10.983/10.912|1515.6/1522.9|-343.446256564846|
|FULL|1764|4/4|30/30|9.623/9.556|1525.7/1530.1|-343.446302545201|

Maximum final-energy difference is2.96e-12 Eh; iterations/HVP counts are
unchanged. Whole-task times and RSS are essentially unchanged in these single
runs. The kernel improvement is **not** a measured complete-SCF speedup here.
The broad source support of these subspaces does not give the new factor
adjoint an asymptotic advantage, and other allocations dominate RSS.

## Remaining work

Accepted orthogonal states/residuals, generator adjoints, directional and
coupled-response product arrays, selected-state producers and restricted-metric
carrier metadata remain. These must become factor/tile consumers before the
complete gradient/HVP pipeline can claim bounded product memory. In particular,
the factor integral-adjoint speedup does not establish a factorized HVP.

Next: differentiate the shared factor graph and implement the unsymmetrized
one-body/factor pullbacks, preserving the residual and state-average gauge
terms. Remove obsolete full-carrier consumers only after their replacements
pass complete gradient/HVP and SCF comparisons. Root plan:
`TODO_ORTHOGONAL_CI.md`; theory: `article/source_column_hamiltonian.tex`.

CSV generation: `tests/mcscf/summarize_factor_gradient.py`. Raw logs retain
input/binary/source hashes, configuration, affinity, timings and failed candidates.
