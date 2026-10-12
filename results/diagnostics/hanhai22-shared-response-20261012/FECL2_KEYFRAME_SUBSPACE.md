# FeCl2: keyframe corrections outside the retained AH space

## New evidence

Slurm 1249453 on Hanhai22 anode31, original FeCl2 guess, 32 bound cores,
exact integrals, unchanged CIAH controls. Five unit executables passed.
The new diagnostic uses cached AH metric images; it evaluates no extra HVP.
The native calculation converged in 14 macros / 281 HVPs, energy
-2181.624464869070 Eh, GNU wall 9.94 s, peak RSS 493.98 MiB.
This is attribution, not a new convergence improvement. No iron long run
was submitted and no threshold or candidate policy was modified.

For trial basis Q and complete physical metric G, the portion of a covector
u captured by the space has squared norm

$$
\|u\|_{Q,G}^{2}=(Q^Tu)^T(Q^TGQ)^{-1}(Q^Tu).
$$

Here u is the true keyframe gradient minus the preceding model gradient,
in the same immutable origin chart. The complete dual norm of u was already
computed. This defines a metric projection, not a Euclidean basis fraction.
The additional storage is only the existing at-most-30-dimensional projected
metric. Independent unit tests verify a non-diagonal metric projection and
that no Hessian image is requested by the diagnostic.

Let s denote displacement since the preceding keyframe. The signed cosine
is u^T s divided by the product of its complete dual norm and s's physical
primal norm. The table gives each keyframe, not just a macro average.

| Macro | Captured correction squared | Signed cosine |
|---|---:|---:|
| 10 | 20.46%, 24.01% | 0.0020, 0.0038 |
| 11 | 0.26% | 0.0017 |
| 12 | 12.63%, 13.99% | 0.0008, 0.0016 |
| 13 | 2.58%, 22.30% | 0.0021, 0.0016 |
| 14 | 9.94%, 13.66% | 0.0005, 0.0008 |

Thus most keyframe correction is outside the sampled trial space. It is
also almost perpendicular to the actual displacement. This is consistent
with slow expansion after model refresh, but does not by itself prove an
operator defect or that adding this correction direction will improve wall
time. The earlier failed RHS-refresh experiments remain relevant evidence.

## Consequence for a possible local curvature update

If the frozen model is H and the keyframe gradient difference is y,
u=y-Hs. A symmetric rank-one correction u u^T/(u^T s) satisfies the secant
relation but divides by the nearly orthogonal contraction observed above.
These measurements do not justify adopting it.

A bounded-norm alternative in physical whitened coordinates is the
symmetric minimum-change rank-two update. In original coordinates, write
w=Gs and d=s^TGs. Its correction is

$$
\Delta H=\frac{u w^T+w u^T}{d}
 -\frac{u^Ts}{d^2}w w^T.
$$

It is symmetric and directly satisfies Delta H times s equals u. In
whitened coordinates it is the minimum-Frobenius-norm symmetric secant
correction. Unlike rank one, it does not divide by u^T s. Its application
requires vector contractions, not another orbital HVP or dense Hessian.

This is a proposed next controlled experiment, NOT implemented production
behavior and NOT a demonstrated speedup. It must consistently update both
retained AH images and sampled coarse-preconditioner images, while future
operator images include the same correction. Updating only the gradient or
only one cache would produce incompatible quadratic models.
Store updates only inside the current macro; no cross-point secant history
is implied. Physical HVPs, coordinates, retraction, tolerances and true-energy
acceptance remain unchanged. The resulting micro model would be an updated
approximation, not the unmodified exact accepted-point Hessian.

Raw evidence: fecl2-ritz-cost-1249453/native/run.err, neo.tsv, output.out,
resource.time and hashes.txt. summarize_ritz_cost.py extracts all keyframes.

## Implemented rank-two experiment and ablation: rejected

The proposed update was implemented and tested, with the same update applied
to cached AH images and to the coarse-preconditioner images. Future coupled
images included both orbital and structure parts of the correction. A
keyframe invalidated the old pending residual and eigenvalue comparison, so
the next AH solve reprojected in the retained space before expanding it.

Independent tests checked symmetry, the exact secant equation for an
orthogonal mismatch/displacement (where SR1 is undefined), agreement of
split and full images, cached AH images against a dense updated operator,
and exact coarse inverse action on updated samples. Five unit executables
passed. No convergence threshold, component bound, microstep budget,
physical HVP, retraction or true-energy acceptance condition changed.

| FeCl2 experiment | Job | Macro / HVP | AH Ritz solves | Energy / Eh | Wall / s | RSS / MiB |
|---|---:|---:|---:|---:|---:|---:|
| Archived retained diagnostic | 1249453 | 14 / 281 | 281 | -2181.624464869070 | 9.94 | 494.0 |
| Rank-two update with reproject | 1249463 | 14 / 261 | 283 | -2181.624464860028 | 9.67 | 524.1 |
| Reproject only, unchanged curvature | 1249468 | 15 / 272 | 295 | -2181.624464861139 | 11.02 | 556.8 |

All three native FeCl2 observations ran on anode31 at 32 bound cores. The
retained run was reused, not rerun. The trial endpoints agree with the
PySCF coupled-Newton reference within about 2.4e-8 Eh. Nevertheless, the
rank-two variant did not reduce outer steps or resolve the micro-budget
tail. Its 20 fewer HVPs coexist with 22 retained-space Ritz evaluations
that needed no additional HVP; action-count differences must not be
interpreted as a 20-direction improvement in model convergence. The
reprojection-only ablation is worse. Single-run wall differences do not
establish a small statistical speedup, and no allocation capture establishes
the cause of the RSS differences.

Frozen rank-two binary regression, array 1249464:

| Case | Retained macro/HVP | Trial macro/HVP | Trial energy / Eh | GNU wall / s | RSS / MiB |
|---|---:|---:|---:|---:|---:|
| F2 | 3 / 23 | 3 / 23 | -198.751155830526 | 0.33 | 17.2 |
| 241 | 3 / 21 | 3 / 19 | -230.720590392726 | 2.75 | 801.9 |
| MnF2 | 10 / 144 | 10 / 135 | -1348.893353224001 | 4.94 | 553.1 |

All converge with unchanged macro counts and energy differences below
1.2e-10 Eh from the retained references. These are separate-node small
regressions, not matched performance claims.

Conclusion: local symmetric secant correction is algebraically consistent,
but this implementation does not demonstrate a material improvement of the
requested convergence efficiency. Its production class, cache-update APIs,
test cases and reprojection-only ablation were removed locally and from
the cluster source mirror. No switch, fallback or cross-point history was
retained. The mirror main executable was restored from frozen 1249453;
restoration did not execute an old benchmark. The original projection
diagnostic and its independent unit test remain. No iron long run was
submitted on the basis of this trial.
