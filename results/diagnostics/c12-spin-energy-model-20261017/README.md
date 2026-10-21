# Selected-space CI and AH spin-energy preconditioning

2026-10-17. C12/LOFLEA HAO CAS(12,12), exact Libcint/cc-pVDZ,
structure-space Davidson, ISCF=9/CIAH. New runs use Hanhai22/anode65,
32 physical cores, cores/close binding and one BLAS thread. Saved controls
used anode42; both nodes have AMD EPYC 9654 CPUs. Timings are single-run,
cross-node observations, not same-node speedup measurements. Old versions
were not rerun. Unrelated local RI edits were not transferred.

## What caused the extra actions

The retained FULL solver uses same-dimensional orthonormal CI coordinates
and a carrier energy-diagonal preconditioner. It does not introduce
independent spin-product coefficients. The broad selected-space solver
instead used an occupation-scaled spin Gram inverse, approximating the
metric but ignoring the energy denominator. Its structure AH correction
still used raw scalar Jacobi denominators, discarding spin recoupling that
was already available in the accepted point. These are weak inverse models,
not changes to the physical selected space or evidence of an incorrect HVP.

In the saved ION first-macro trace, the structure residual contributes
99.7% of the squared coupled residual (orbital 0.04669, structure 0.84907).
The saved HVP counts are 12,12,14,18,27,3, versus FULL's 8,7,22,3.
Thus the extra HVPs are not explained solely by terminal overprecision.
The final warm CI counts, 22 versus 6, are not cold-start iteration counts.

For arbitrary selected columns W, G=W^T W is block diagonal by occupation.
The new shared inverse model is P=K G K, K_ii=sqrt(d_i/G_ii), with positive
d matching the existing energy-diagonal model entry by entry. Its inverse
uses the retained G factors and two diagonal scalings. CI uses
d_i=abs(h_i-E m_i); AH uses abs(h_i-(E+lambda)m_i). Both retain their
existing numerical denominator floors. The signed physical Hessian and
negative-curvature treatment are unchanged.

For the AH tangent, A=M C contains the existing root constraints. The
restricted inverse is R-R A(A^T R A)^-1 A^T R, where R=P^-1, followed
by the existing tangent-coordinate projection. This follows by eliminating
the multiplier in P u+A xi=r, A^T u=0. Merely projecting R is also a valid
positive preconditioner, but is not the inverse of the restricted model.
The selected-space physical metric still requires its original solve;
a selected block of the FULL inverse is not substituted for it.

Only O(n_str*n_state) temporary vectors and a small n_state-square factor
are added. Existing occupation factors are shared, not rebuilt or copied
per thread. No dense structure Hessian, global spin-pair cache, or extra
carrier image is introduced. No model branch names ION, C12 or this input.
Physical actions, orbital coordinates, CI transport and certification targets
are unchanged. The transport recycling experiment remains removed.

Derivation: article/source_column_hamiltonian.tex, Section 10.9,
"A diagonal-matched spin model on the CI tangent". Its PDF was rebuilt
and the added equations inspected.

## Independent checks and same-operator diagnosis

The new unit checks compare the inverse model with a separately constructed
small dense Gram, and compare the actual AH implementation with the inverse
of the independently restricted model for one and two selected states.
These use a nonclosed selected space. Existing physical action, gradient,
HVP and transport tests are unchanged. Qualification 1265250 passed 11/11
for AH alone; qualification 1265263 passed 11/11 for CI+AH (3.50 s).

Diagnostic 1265254 uses the saved endpoint orbitals but the current operator,
with identical cold-start seeds and stopping targets for each inverse model.
It is not an old-code SCF rerun or a repeated fixed-point energy-validation
exercise. The reported CI energy is active-space energy, not total energy.

| Current-operator cold CI model | H/S block actions | CI residual / Eh | Wall / s |
|---|---:|---:|---:|
| Metric-only spin inverse | 100 | 8.8275e-8 | 8.9723 |
| Discarded occupation-trace model (diagnostic only) | 79 | 9.7924e-8 | 6.9067 |
| Diagonal-matched positive spin model | 81 | 9.7058e-8 | 7.0521 |

All three active energies agree within 1.5e-14 Eh. At this point no energy
diagonal or occupation-trace denominator is negative. Their range is
10.0644--38.4174 Eh; maximum within-block normalized spread is 1.2024.
Therefore the earlier combined CI/transport experiment cannot be dismissed
as simply a negative-denominator failure. It changed multiple ingredients
and did not isolate their causes. The discarded trace formula is retained
only in this diagnostic, never as a production switch or fallback.

For one deterministic structure-Hessian RHS, the best scalar linear residual
after preconditioning falls from 0.568309 (raw Jacobi) to 0.397262 (the new
constrained spin model). This is an inverse-model diagnostic, not a claim
of a proven condition-number bound or an exact Hessian inverse.

## End-to-end convergence

The ION inputs are byte-identical, SHA256
48b44d64b773d4abaf388fcf728d478a09ab76f5426f889822b1dbc8e031432f.
Original guesses and ITMAX=2000 are unchanged. Energy tolerance remains
1e-7 Eh, projected gradient tolerance 1e-3, CI stationarity 1e-7 Eh and
metric-CG true relative residual sqrt(machine epsilon).

The saved FULL result is a different variational space, not an energy-parity
target for ION. Timing columns distinguish application end-to-end time from
GNU-time process wall time. Peak RSS is the process maximum, converted from
KiB to GiB; it is not Slurm's intermittently sampled RSS.

| Version / job | Space | Structures | Steps | HVP | Final CI actions | Final total energy / Eh | Application / s | Process / s | Peak RSS / GiB |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Saved FULL / 1265114_2 | FULL | 226512 | 4 | 40 | 6 | -399.417495587055 | 102.727418 | 112.47 | 3.572872 |
| Saved selected metric / 1265216_1 | ION(0:2) | 47256 | 6 | 86 | 22 | -399.397232509246 | 129.879149 | 132.79 | 2.932171 |
| New AH only / 1265253_1 | ION(0:2) | 47256 | 4 | 44 | 24 | -399.397232509248 | 111.176418 | 114.47 | 2.947399 |
| **New CI+AH / 1265266_1** | **ION(0:2)** | **47256** | **4** | **44** | **19** | **-399.397232509245** | **104.286136** | **107.11** | **2.823959** |

Against the saved ION baseline, observed application time is 19.71% lower,
HVP count is 48.84% lower, and peak RSS is 110.81 MiB lower. The final
energy differs by +1.02e-12 Eh. These are completed, converged runs, not
one-step projections. The final ION input hash and CPU affinity 96-127
match the saved ION run, but the nodes differ as stated above.

The two new runs are on the same node with the same affinity. Adding the
CI inverse model to the constrained AH change leaves 4 steps and HVP
counts 9,5,27,3, reduces final warm CI actions from 24 to 19, and reduces
observed application time by 6.20%. This distinguishes the large HVP benefit
of the AH change from the smaller CI benefit. The full CI+AH endpoint has
projected |g|_inf=3.16706342e-7 and CI residual 9.27424325e-8 Eh.

This does not solve every CI conditioning issue: final ION CI still needs
19 actions versus FULL's 6. The model keeps occupation-local spin coupling
and energy scales but does not represent all orbital-induced cross-occupation
couplings of the restricted problem. The data do not establish that this
remaining omission alone accounts for all thirteen actions. No new FULL
run or new threshold was used to obtain the comparison.

The final frozen bundle is c12-spaces-bundle-1265264, qualified by
1265263. Its executable SHA256 is
d6dde8f05ac10270dff15e079febcfd4fc8d21b94bd078fb9662af8a3cc3ec3c.
The AH-only bundle is c12-spaces-bundle-1265252, SHA256
9362c60e786a35408837c713bb55b6532b2c6a8ca29f926681bc75b18cdc2e54.
Source manifests were checked before freezing, and retained production
source hashes match the final qualification manifest. All new native jobs
completed with exit 0. No long optimization remains running from this round.

Final numerical records are in comparison.csv; full native outputs, inputs,
traces, manifests and resource logs are in evidence.tar.gz. Executables
remain in the frozen Hanhai22 bundles and are not copied into Git.
