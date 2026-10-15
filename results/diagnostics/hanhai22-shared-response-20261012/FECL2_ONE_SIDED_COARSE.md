# FeCl2: one-sided coarse correction versus balanced enrichment

Rejected experiment; the current source retains the symmetric balanced inverse.

## Motivation

Native 1249619 localizes the tail residual almost entirely to the orbital
sector. Keyframe corrections are mostly outside the actual AH trial space.
The sampled coarse space is larger: every coupled direction supplies
separate orbital and structure directions/images without another HVP.
Its balancing projector acts on this larger space, not just the AH space.
Test whether that second projection suppresses useful enrichment before
changing any stopping criterion or adding curvature evaluations.

## Algebra

Let L=H-omega G denote the shifted coupled operator, Q the actual AH basis,
V the sampled split basis, E=V^T L V, and P the unchanged base inverse.
Write K=VE^{-1}V^T. The retained balanced inverse is

$$
P_{\mathrm{bal}}=K+(I-KL)P(I-LK).
$$

The diagnostic experiment instead uses

$$
P_{\mathrm{one}}=K+P(I-LK).
$$

Both obey the exact sampled inverse condition

$$
P_{\mathrm{bal}}LV=V,\qquad P_{\mathrm{one}}LV=V.
$$

Their difference maps into the sampled coarse span:

$$
(P_{\mathrm{one}}-P_{\mathrm{bal}})r
=KLP(I-LK)r\ \in\operatorname{span}V.
$$

If Q spans V, subsequent Davidson orthogonalization eliminates this
difference. Here Q is a proper subspace of V, so split-sector information
in the remaining coarse directions can change enrichment. This is not the
previous experiment that removed orbital--structure coupling from E, or
the experiment that doubled every actual Ritz direction. E and the AH
workspace remain unchanged.

P_one is generally nonsymmetric. It is used only to generate a trial
direction; the actual Hessian/metric and their generalized Ritz projection
remain symmetric and certified. Do not claim SPD-PCG convergence bounds
for this indefinite, nonlinear CIAH experiment. General background on
balancing and deflation is provided by
[Tang et al., Journal of Scientific Computing (2009)](https://link.springer.com/article/10.1007/s10915-009-9272-6).
The above identities are verified independently rather than relying on a
literature theorem outside its assumptions.

## Cost and controls

For each residual, solve Ea=V^Tr and set

$$
d=P(r-LVa)+Va.
$$

This removes one small coarse solve and one tall projection from the
balanced formula. It still makes one base inverse call, no additional
HVP/metric action, and allocates no additional large workspace.
HVP, nonorthogonal coordinates, retraction, phase handling, analytic
diagonal, micro-budget, keyframes and convergence thresholds remain fixed.
No new numerical parameter, production selector or fallback is introduced.

Dense tests cover positive/indefinite H, nonidentity G, eigenvalue shifts,
coarse dependence, exact sampled inverse, one base action, and equality of
the complementary projected directions when Q spans V. The existing full
AH/coordinate/diagonal tests remain part of the compute-node gate.

## Completed experiment: rejected

Hanhai22 1249633 completed on anode31; all five units passed, followed by
normal FeCl2 convergence. The retained reference is archived, not rerun.

| Variant / job | Macros | HVPs | Energy / Eh | GNU wall / s | GNU RSS / MiB |
|---|---:|---:|---:|---:|---:|
| Balanced, 1249596 | 12 | 251 | -2181.624464871899 | 8.30 | 516.23 |
| One-sided, 1249633 | 13 | 241 | -2181.624464850508 | 8.47 | 528.17 |

The new energy is +3.3225e-8 Eh above the independent PySCF coupled
Newton endpoint. HVP count is slightly lower, but outer count increases
and total wall time is not improved. All final four macros still end at
the micro-budget. This is not the requested robust tail improvement.
The experimental implementation and its changed tests were removed;
the balanced formula and its symmetry test are restored. No switch or
fallback is retained. The result establishes that the second coarse
projection affects enrichment but is not its dominant cause.
