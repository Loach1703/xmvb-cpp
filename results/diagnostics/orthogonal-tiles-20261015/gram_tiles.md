# Selected Gram transport: structural tile rewrite

## Implemented scope

The orthogonal-CI branch of RepresentationOverlap now evaluates the prepared
metric, cross overlap, directional overlap, orbital pullback, and pullback
direction through selected-product panels. The old full-carrier helpers in
representation.cpp were deleted, not retained as an alternative implementation.
The physical non-complete-string pair representation is unchanged.

The beta pullback uses the input-side exterior generator identity. It does not
require a full product-space transition density or neighboring output panels.
Only original touched beta columns enter its seed contraction. Alpha and beta
are handled independently, including unequal electron counts. See
`article/selected_gram_tiles.md` for the derivation and storage bounds.

This does not migrate the separate accepted-state integral adjoint caches or
the full structure-response sigma/adjoint interface. It does not claim that
all VBSCF consumers are now free of full product arrays.

## Final production regression

Job 1260142: Hanhai22, anode15, 32 cores, exact libcint, cc-pVDZ, NEO,
Davidson, analytic Hessian diagonal. The original input orbitals and thresholds
(gradient 1e-3, energy 1e-7 Eh) are unchanged. All five spaces ran in one Slurm
allocation. Six unit regressions passed before native runs.

Historical job 1260033 ran on anode01. That node was requested for the new run,
but was fully allocated with a scheduler estimate two days away; the pending
job was moved to anode15. Therefore these are cross-node historical timings,
not a matched-allocation speedup measurement. No historical executable was run.

| Space | Structures | Steps old/new | HVP old/new | Old SCF s | New SCF s | Old process s | New process s | Old RSS MiB | New RSS MiB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| COV | 14 | 3/3 | 20/20 | 2.679260 | 2.118234 | 4.28 | 3.68 | 1572.56 | 1537.23 |
| ION(0:1) | 294 | 3/3 | 27/27 | 3.765359 | 3.537664 | 5.32 | 5.11 | 1566.61 | 1546.43 |
| ION(0:2) | 1134 | 4/4 | 32/32 | 4.515104 | 4.470734 | 6.03 | 6.00 | 1606.10 | 1548.15 |
| ION(0:3) | 1694 | 4/4 | 32/32 | 4.691147 | 4.632540 | 6.34 | 6.23 | 1622.50 | 1543.91 |
| FULL | 1764 | 4/4 | 30/30 | 3.545873 | 3.387353 | 5.15 | 4.90 | 1572.07 | 1566.88 |

| Space | New final energy / Eh |
|---|---:|
| COV | -343.239211501261 |
| ION(0:1) | -343.407736636062 |
| ION(0:2) | -343.443415326851 |
| ION(0:3) | -343.446256564846 |
| FULL | -343.446302545201 |

Maximum absolute final-energy difference from job 1260033 is approximately
2e-12 Eh. Native SCF steps and HVP counts are unchanged. RSS remains dominated
by the molecular calculation outside these small carrier workspaces; this
test does not demonstrate a large end-to-end memory reduction.

## CAS(12,12) transport kernel validation

Job 1260182: same final binary sources, anode15, 32 cores. These are synthetic
fixed-orbital transport kernels, not molecular energies, H/S sigma timings,
or converged SCF calculations. There is no saved matched old transport-kernel
timing, so no kernel speedup ratio is reported.

| Structures | Touched original products | Reachable products | Gram s | Tangent s | Pullback s | Pullback tangent s | Peak RSS MiB |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 64 | 136374 | 0.10620 | 0.12890 | 0.09674 | 0.21840 | 37.63 |
| 3 | 128 | 157108 | 0.11283 | 0.07350 | 0.07020 | 0.19652 | 39.39 |
| 20 | 376 | 192380 | 0.13138 | 0.07506 | 0.07843 | 0.29537 | 37.89 |
| 132 | 924 | 226512 | 0.13593 | 0.15970 | 0.12668 | 0.31654 | 38.42 |

Largest action/adjoint duality error: 4.03e-16. Largest relative pullback
direction finite-difference error: 1.79e-8, with step 1e-6 and acceptance 2e-7.
Independent occupied-minor unit tests additionally cover sparse arbitrary
nonclosed spaces, unequal spins, identity charts, and narrow/wide panels.

## Interpretation

The change removes mandatory full product images from this transport chain
and preserves its derivatives. It is not an order-of-magnitude SCF speedup.
In particular, ION(0:3) remains slower than FULL. Restricted metric solves,
other response consumers, compound-row construction, and structurally broad
transformed support remain relevant. One CAS(12,12) structure already reaches
136374 products, versus 226512 for all 132 covalent structures; no factor of
132 in runtime follows from the structure count alone.

Raw evidence: `240-metric-1260142/`, `gram-tiles-1260182/`.
Machine-readable production comparison: `240_gram_results.csv`.
Intermediate validation jobs 1260100, 1260121 and 1260139 tested the staged
implementation; they are not used as historical production baselines.
