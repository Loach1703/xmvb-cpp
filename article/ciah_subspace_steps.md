# Separating CIAH subspace enrichment from microstep advancement

## Scope

The orbital chart, finite CI transport, analytic HVP, generalized metric,
preconditioner and outer stopping conditions are unchanged. This note concerns
the coupled orbital--structure vector in one immutable accepted-point chart.
No closure of an arbitrarily selected VB structure space is assumed.

Let $H_0$ and $G_0$ denote the matrix-free Hessian and metric actions. At a
keyframe $p_a$, the true pulled-back model gradient is $g_a$. The current
displacement from that keyframe is $s=p-p_a$, with predicted gradient

$$
g=g_a+H_0s.
$$

The retained basis $Q$ has cached images $Y=H_0Q$ and $Z=G_0Q$, and small
projections $B=Q^TY$, $M_Q=Q^TZ$. Enrichment adds one column and one HVP.
Changing the gradient only changes $Q^Tg$; it does not require new images.

## What the PySCF reference does

PySCF's `soscf.ciah.davidson_cc` reads a gradient callback, solves the small
augmented problem and yields a candidate before expanding the Davidson space.
`mcscf.newton_casscf.update_orb_ci` independently decides whether to consume
that candidate. Thus a nonconverged Ritz step is not inherently an error.
After a nonconverged yield is resumed, the generator normally enriches using
the yielded residual; a converged yield may resume without another HVP.

References inspected on 2026-10-17:

- https://pyscf.org/_modules/pyscf/soscf/ciah.html (`davidson_cc`)
- https://pyscf.org/_modules/pyscf/mcscf/newton_casscf.html (`update_orb_ci`)

We do not transplant the reference implementation's absolute start tolerance
or its orthogonal orbital representation into the nonorthogonal VB chart.
Our existing cached reprojection already avoids expanding an obsolete RHS.
However, its consumption rule requires relative forcing for a cached-only
candidate while allowing a newly enriched candidate without that forcing.

## Why an accurate AH increment need not be an efficient total step

Write the current local displacement as $s=Qz$. An exactly solved projected
AH increment $\delta=Qd$ with root $\omega$ satisfies

$$
(B-\omega M_Q)d=-Q^Tg_a-Bz.
$$

The cumulative coefficients $z+d$ therefore obey

$$
(B-\omega M_Q)(z+d)=-Q^Tg_a-\omega M_Qz.
$$

The last term depends on previous increments. Tightening the Ritz residual
does not remove it. Conversely, repeatedly replacing the trajectory by one
fixed-anchor AH target can stall at a damped point: for $H=G=g_a=1$, the AH
step is $(1-\sqrt{5})/2$, whereas the Newton step is $-1$. Neither blindly
tightening the solve nor freezing a single damped AH target is the proposal.

## Positive interior correction from cached images

If $B$ is positive definite, the quadratic restricted to the current span has
the unique correction

$$
d_N=-B^{-1}Q^Tg,\qquad \delta_N=Qd_N,
\qquad H_0\delta_N=Yd_N,\qquad G_0\delta_N=Zd_N.
$$

These use a small Cholesky factor reused until the basis changes. No HVP,
CI response, physical metric inverse or full-space array is added. The metric
does not enter the unshifted Newton equation; it remains in the AH equation,
step measurement and the existing physical geometry.

For a nonzero correction,

$$
Q^T(g+H_0\delta_N)=0,\qquad
g^T\delta_N=-\delta_N^TH_0\delta_N<0.
$$

If $s\in\operatorname{span}(Q)$, the new total local point is

$$
s+\delta_N=-QB^{-1}Q^Tg_a.
$$

Thus the same-span total target is independent of earlier damped increments.
At the same time, the complete residual

$$
r_N=g+H_0\delta_N
$$

is retained. Its components outside the span need not be small; projected
stationarity is neither full Newton convergence nor nonlinear convergence.
In particular, the full gradient norm can grow despite quadratic descent.

## Consumption and enrichment

1. Form the usual AH candidate with the existing cached reprojection/enrichment
   process. Keep the existing minimum start dimension.
2. If the positive projected correction has nonzero predicted decrease and
   satisfies the *existing* orbital and CI component cap without clipping,
   consume it even if the complete relative AH forcing test is unresolved.
3. Consume that Newton correction only once per pair of basis dimension and
   accepted true-keyframe generation. A second correction in the same span
   and same linear model would only remove floating-point residue.
4. A changed predicted RHS is reprojected before any further enrichment. If
   the current span cannot provide another usable correction, its complete
   residual generates the next Davidson direction using the existing
   preconditioner. Do not certify only the projected residual.
5. Nonpositive projected curvature and steps outside the existing component
   cap remain AH cases. This preserves negative-curvature/globalization
   behavior instead of inverting an indefinite quadratic.
6. Existing gradient-growth rollback, keyframe update, microstep budget,
   finite orbital/CI update and relaxed physical-energy acceptance remain.

The generation counter is model state, independent of diagnostic evaluations.
The Newton branch is recorded explicitly: root zero and reference component
one describe its unshifted residual, not a converged AH eigenpair.

This is a direction-selection change in the positive/interior regime, not
an algebraically identical rearrangement of all CIAH trajectories. The
factorization is not a condition-number certificate; nonfinite corrections
are unavailable and infeasible corrections remain in the existing AH regime.

## Validation and interpretation

Unit checks cover cached action counts, a nonidentity metric, exact projected
stationarity with a nonzero complement, indefinite/singular curvature, zero
corrections, and current-RHS reprojection before enrichment. Native production
validation uses the original iron-porphine singlet input on Hanhai25, 32
physical cores, exact Libcint/cc-pVDZ, FULL OEO CAS(8e,11o), analytic diagonal,
gradient $10^{-3}$ and energy-change $10^{-7}$ Eh. No old calculation is rerun.

Historical 17/242 and current 21/316 trajectories do not stop at exactly the
same point: their energies differ by about $4.53\times10^{-6}$ Eh. Existing
22/318 data reach the same lower endpoint as 21/316. Report both comparisons,
including true energies, macrosteps, HVPs, wall time and peak RSS. Fewer small
Ritz solves alone do not demonstrate a speedup or better convergence.
