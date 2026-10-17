# Advancing a clipped AH proposal without additional Hessian actions

## Scope

This is a bounded experimental change to the CIAH microstep driver, not a
change to the orbital chart, finite orbital/CI transport, Hessian action,
metric, preconditioner, root selection or convergence thresholds. Validation
starts with small systems; the previous split-image enrichment experiment
is not reintroduced.

The path only scales an already admitted coupled direction. It does not
assume complete-CAS closure or replace the selected VB structure space.
Measured results and provenance are recorded in
`benchmarks/ciah-path-20261017/README.md`.

## Proposal and path

Let $g_a$ be the full gradient at the point where an augmented-Hessian (AH)
candidate $t$ is generated. In the retained trial space $Q$, the candidate
obeys the projected generalized equation

$$
Q^T\bigl[g_a+Ht-\omega Gt\bigr]=0,
$$

where $H$ and $G$ are the frozen Hessian and metric actions and $\omega$ is
the selected AH root. The workspace also retains the complete residual,
including components outside $Q$; projected stationarity is not full
stationarity.

If the component limit clips $t$, the original driver changes its predicted
gradient and obtains another AH candidate before completing the first one.
That candidate can require another subspace expansion. The proposed change
instead treats the admitted candidate as one finite path. With the remaining
fraction $\rho_j$, its next increment is

$$
\delta_j=\alpha_j t,\qquad
\alpha_j=\min\!\left(\rho_j,\frac{d_{\max}}{m(t)}\right),\qquad
\rho_{j+1}=\rho_j-\alpha_j,
$$

where $m(t)$ is the existing maximum orbital/structure component measure and
$d_{\max}=0.03$ is the unchanged component bound. The measure is homogeneous,
so $m(\alpha_jt)=\alpha_jm(t)$. The matched image is already known:

$$
H\delta_j=\alpha_jHt.
$$

No additional Hessian action, metric action, preconditioning or Ritz solve is
required to form a path increment. If the entire path is completed,

$$
\sum_j\delta_j=t,\qquad \sum_jH\delta_j=Ht.
$$

Even if the path is interrupted, its partial sums obey exactly the same
frozen model relation, up to arithmetic roundoff:

$$
g_{j+1}^{\mathrm{pred}}=g_j^{\mathrm{pred}}+H\delta_j.
$$

## What the path does not certify

The source AH eigenpair does not become a new AH eigenpair for the gradient
at an intermediate point. For example, if $g_j=g_a+(1-\rho_j)Ht$, the residual
of the remaining direction at the old shift is

$$
g_j+H(\rho_jt)-\omega G(\rho_jt)
=g_a+Ht-\omega Gt+\omega(1-\rho_j)Gt.
$$

The last term is generally nonzero. Consequently, replayed pieces cannot
inherit an AH convergence/forcing flag for the current right-hand side.
The recorded source Ritz residual remains diagnostic only; the complete
predicted gradient still controls the existing microstep/keyframe tests.

## Invalidation and enrichment

Every piece still executes the original gradient-growth, microstep-budget,
keyframe, and subspace-exhaustion checks. A kept true keyframe clears the
unused path immediately, since its corrected gradient can change the useful
direction. No additional keyframe rule is introduced.

When the path is finished or invalidated, the workspace first reprojects the
current complete gradient using its cached images. Subsequent expansion
uses the residual of that current equation. An outside-space residual is
therefore delayed while an existing proposal is being followed, not dropped
or declared converged. Physical macro acceptance and final energy/gradient
tests are unchanged. The existing exhaustion stop can terminate a clipped
proposal before completion; this experiment does not expand that budget.

The path borrows the driver's existing proposal and its operator image;
only its remaining fraction is new state. The proposal is not replaced
while a path is active. No per-thread AO-integral or spin-product arrays
are allocated. The trace
distinguishes `replayed`, `new_hvp`, and `residual_rhs_current`. The existing
`micro_iterations` counter counts Ritz solves; the `ciah_micro` events count
actual consumed pieces, including replayed ones.

The proposal targets action count in clipped phases. It does not promise
improvement in an unclipped tail or assume that fewer HVPs imply faster
end-to-end convergence. Changed nonlinear trajectories must be measured.
