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
