# FeCl2: AH expansion and exact active-source row blocks

## What the retained measurements establish

The native scalar diagnostic 1249619 (Hanhai22 anode04, 32 bound threads)
converged in 12 macros / 249 HVPs, at -2181.624464872068 Eh, in 7.40 s
with peak RSS 523.69 MiB. This is archived evidence, not a rerun.

| Instrumented scope | Seconds |
|---|---:|
| AH augmented spectrum and image reconstruction | 0.01977 |
| AH projection update | 0.02421 |
| Complete orbital dual norms | 0.23973 |
| Complete structure dual norms | 0.08921 |
| HVP expansion, including observation | 1.93116 |
| Analytic diagonal builds | 0.82594 |

These are nested scopes, not an exhaustive wall-time partition. Do not
add the complete norm scopes to Ritz residual norm time. The structure
Riesz solves require at most one CG iteration and have certified relative
residual below 8.85e-16. The maximum relative Galerkin residual is 5.11e-12.
Neither a loose structure norm nor an inaccurate small projected solve
explains the tail. The last budget-exiting macro has 99.91% of the residual
squared in the orbital sector. Keyframe mismatches in that macro have only
3.43% and 3.17% of their squared magnitude captured by the current AH basis.
This is evidence about sampled directions, not a full Hessian spectral bound.

PySCF's coupled CASSCF implementation accepts usable approximate AH
microsteps, updates the predicted gradient, and checks its gradient-dependent
microstep budget before deciding whether to refresh the keyframe. Its
Davidson generator retains trial vectors and Hessian images as the gradient
changes. Thus a budget exit is not by itself a missing implementation, and
raising the budget does not address ineffective complementary directions.
See the pinned [coupled CASSCF implementation](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/mcscf/newton_casscf.py)
and [AH Davidson implementation](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/soscf/ciah.py).

## Rejected virtual-Fock approximations

Both original-guess trials used anode31, 32 bound threads, exact integrals,
the same HVP/chart/retraction and numerical controls. All five unit
executables passed before chemistry. No archived baseline was rerun.

| Variant | Macros | HVPs | Energy / Eh | Wall / s | RSS / MiB |
|---|---:|---:|---:|---:|---:|
| Archived 1249596 | 12 | 251 | -2181.624464871899 | 8.30 | 516.23 |
| Virtual-Fock block 1249720 | 16 | 221 | -2181.624464698809 | 9.37 | 519.20 |
| Inverse-error fit 1249733 | 13 | 261 | -2181.624464860377 | 10.06 | 534.10 |

The first trial couples virtual rows at each occupied source through
2 n_q F_vv, retaining the original analytic diagonal. It is only a partial
kinematic model, not an exact Hessian block. Its energy is 1.84924e-7 Eh
above the strict coupled-Newton PySCF reference -2181.624464883733 Eh.
Fewer actions did not deliver a valid efficiency/precision improvement.

The second trial uses the exact already computed image of an admitted
direction to fit a convex combination of the diagonal and Fock-block
inverses. Let L=H-omega G and let P_D/P_F include the existing symmetric
coarse balancing. For each previously unlearned q, define

$$
e_q=P_D Lq-q,\qquad d_q=(P_F-P_D)Lq.
$$

The current-point least-squares weight is

$$
\tau=\operatorname{clip}_{[0,1]}
\left(-\frac{\sum_q e_q^TGd_q/(q^TGq)}
{\sum_q d_q^TGd_q/(q^TGq)}\right).
$$

It adds no HVP, but costs additional inverses and a metric action. It
does not improve actions, macros or time. Both implementations were
removed before the next trial, with no runtime option retained.

## Exact same-virtual-row active-source blocks

This is a different model: it couples active *source* indices at a common
virtual row. The formula is derived from the fixed-C orbital action in
PySCF's `gen_g_hop`, retaining the one-RDM and contracted two-RDM terms.
It is a preconditioner model in an accepted orthonormal computational frame;
it does not replace the nonorthogonal VB orbital parametrization or update.

Let a,b,u,v,w,z be active indices and p an external virtual index. Let
gamma and Gamma denote the ERI-symmetric spin-free active RDMs used by
the analytic diagonal, and V_c the inactive Coulomb/exchange potential.
Define the active generalized-gradient matrix

$$
\Lambda_{ab}=
\sum_u(h+V_c)_{au}\gamma_{ub}
+\sum_{uwz}(au|wz)\Gamma_{wzub}.
$$

With J_(p,uv)=(pp|uv) and K_(p,uv)=(pu|pv), the fixed-C principal block is

$$
B_p(a,b)=2\left[
(h+V_c)_{pp}\gamma_{ab}
-\tfrac12(\Lambda_{ab}+\Lambda_{ba})
+\sum_{uv}J_{p,uv}\Gamma_{uvab}
+\sum_{uv}K_{p,uv}(\Gamma_{auvb}+\Gamma_{uavb})
\right].
$$

The factor two matches antisymmetric orbital-coordinate packing. For
a=b this reduces exactly to the existing analytic Hessian diagonal:
the two-RDM term is J Gamma + K W, with
W_uv=Gamma_auva+Gamma_uava. The existing density contraction computes
J(Gamma+2W)-G[W], and G[W]=2J[W]-K[W].
No structure-response contribution is inserted into this orbital block;
the actual coupled HVP still contains the full orbital/structure action.
Off-row and core/active-source couplings remain outside this approximation.

## Construction, memory, and verification gates

The same batched symmetric AO operator supplies K without a full
four-index MO transform. For a symmetric active pair density D_uv whose
off-diagonal entries are one half,

$$
\langle Q_p|G[D_{uv}]|Q_p\rangle=2J_{p,uv}-K_{p,uv}.
$$

Therefore all packed active channels recover the needed K slices. The
shared AO channel workspace keeps the existing memory budget. No thread
receives a full integral graph, Hessian, or row-block collection.
The increased density-channel cost must be measured, not hidden.

| Construction/storage | FeCl2 | Iron porphine |
|---|---:|---:|
| Active orbitals | 8 | 11 |
| Virtual frame rows | 103 | 369 |
| Original core + diagonal-density channels | 42 | 112 |
| Core + complete packed-active channels | 61 | 155 |
| Row-block doubles | 6,592 | 44,649 |
| One row-block collection / KiB | 51.50 | 348.82 |

Signed row eigendecompositions are shared once per accepted point;
storage is O(N_virtual n_active^2), factorization O(N_virtual n_active^3),
and each inverse application O(N_virtual n_active^2), plus the unchanged
source/frame transforms. Two spectral arrays replace the small row-block
collection in the chart; temporary construction data is also shared.
No full Hessian is claimed or stored.

Required tests: independent exact/RI four-index slice references; equality
of block diagonals to the old analytic formula; mixed energy finite
differences for two active-to-same-virtual directions; signed dense
pullback inverse agreement in a nonorthogonal full-AO chart. Molecular
adoption requires unchanged reference-energy precision and better total
cost. No iron long run follows a failed FeCl2 gate.

## Native result: reject this preconditioner

Slurm 1249750, anode31, passed all five unit executables before chemistry.
The independent mixed finite-difference error was 2.77047e-9; the
four-index-reference diagonal errors were 7.11e-15 (exact) and 5.55e-16
(RI). The signed nonorthogonal dense inverse/pullback test also passed.

The original FeCl2 guess converged in **29 macros / 324 HVPs**, at
**-2181.624464866728 Eh**, with **13.87 s** GNU wall and **507.07 MiB**
peak RSS. The energy is 1.7005e-8 Eh above the strict PySCF reference,
but the cost and macro count are inferior to the archived baseline.
The final macro still exits at the microstep budget: 26 microsteps,
model-gradient norm 1.79843e-4, zero clipped steps; 99.64% of the last
residual squared is orbital. No budget/threshold was changed.

| Candidate scope | Seconds |
|---|---:|
| HVP expansion, including observation | 2.77918 |
| Analytic model builds | 2.60279 |
| Complete orbital dual norms | 0.39581 |
| Complete structure dual norms | 0.13730 |
| Small AH spectrum/image reconstruction | 0.02195 |

The largest relative Galerkin residual is 6.97e-12. These results again
exclude the small eigensolve and structure metric solve as the explanation
for stalled AH expansion. Correct principal blocks do not certify a good
complementary-space inverse for the full coupled operator: omitted core,
off-row and structure couplings still matter. This experiment does not
prove which one dominates, and no molecule-specific modification is used
to rescue it.

The trial production fields, kernel, row eigenspaces and tests were removed
as a whole. No optional/fallback preconditioner remains. The original
retained operator and controls are restored; there is no claimed speedup.
Only this derivation/evidence and a per-macro parser for existing native
timing events are retained. No archived molecular baseline was rerun and
no iron long run was submitted.

Restored-source build-only job **1249759** completed on anode31. All five
unit executables passed with scalar tracing. This gate did not execute
another molecular baseline. The parser's pre-existing aggregate fields
were checked to remain identical to HEAD on four archived/native logs;
the added per-macro rows reproduce their accepted macro/HVP totals.
