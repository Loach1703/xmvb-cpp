# Separating CIAH subspace enrichment from microstep advancement

This note records derivations and unsuccessful performance experiments, not
an active production branch. Both the positive-interior and split-image
prototypes passed algebra tests but failed the iron-porphine end-to-end
target. Their production code was removed. The default solver, physical
actions and stopping thresholds remain those preceding these experiments.

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

## Reuse of the separately computed coupled images

For a physical direction $x=(p,q)$ the current joint action already evaluates

$$
x_o=\begin{pmatrix}p\\0\end{pmatrix},\qquad
x_c=\begin{pmatrix}0\\q\end{pmatrix},\qquad x=x_o+x_c,
$$

$$
H_0x_o=\begin{pmatrix}\mathcal A p\\\mathcal Bp\end{pmatrix},\qquad
H_0x_c=\begin{pmatrix}\mathcal B^\dagger q\\\mathcal Cq\end{pmatrix}.
$$

Both images are available before their sum is returned. The existing coarse
preconditioner already consumes these two samples, whereas the AH workspace
has retained only $x$ and $H_0x$. Consequently, the AH representation has
discarded independently known curvature information. The additional selected
VB space is not replaced by a spin-product space; these are two directions
in the *same coupled orbital/selected-structure coordinate space*.

The image-reuse prototype retains $x_o$ and $x_c$, not $x$ followed by
$x_o$. Their span is identical, but the latter representation can subtract
nearly collinear vectors when one block dominates. Its normalized images
then amplify cancellation error. This caused a metric-symmetry rejection in
the F2 regression; the symmetry threshold was not changed.

Starting from block-pure columns, orthogonalizing a physical trial $x$
against the retained basis makes each block separately orthogonal to every
existing column in that block. Cross-block dot products are zero. Retaining
the two pure blocks preserves this property by induction. Both block images
are supplied directly by the physical action; neither is obtained by
subtracting the other from the combined image.

Each direction, Hessian image and metric image uses identical two-pass
orthogonalization coefficients. Therefore no additional physical action is
needed. For a normalized source $v$ and retained orthonormal basis $Q$, set

$$
c=Q^Tv,\qquad w=v-Qc,\qquad
H_0w=H_0v-(H_0Q)c,\qquad G_0w=G_0v-(G_0Q)c.
$$

After both passes, all three arrays are divided by $\|w\|$. Dependence is
tested using the original `lindep` and the Gram eigenvalue ratio

$$
\frac{1-\|c\|}{1+\|c\|}
=\frac{\|w\|^2}{(1+\sqrt{1-\|w\|^2})^2}.
$$

Before normalization, a block whose squared norm is at most `lindep` times
the squared physical-trial norm is excluded by the existing Gram-rank
criterion. At least one of two orthogonal blocks survives this test.
The right-hand expression avoids subtracting nearly equal numbers. The
operator symmetry and projected metric rank checks still include all columns
and both matrix orientations; no assumed symmetry conceals inconsistent
block images. The existing coarse observer remains independent of AH sample
selection. A preconditioner's optional filtering does not define the solver
space.

The action limit remains thirty *physical evaluations*, not thirty cached
columns. At most two columns per action are admitted, so storage is bounded
by sixty coupled-coordinate vectors and their two images. No AO integrals,
structure-pair array or complete spin-product carrier is copied per thread.
The extra small-matrix eigensolve and basis storage are real costs and must be
included in wall time and RSS. The three-action start rule is unchanged.

All later Ritz solves and their complete residuals use the enriched span.
Thus this increases available exact curvature *per HVP*, rather than asking
the original narrower span to satisfy a tighter residual threshold. It does
not guarantee a factor-of-two HVP reduction: dependence, the shifted root,
clipping, keyframes and nonlinear globalization still matter. The unchanged
reference-weight filter also means that merely exposing a pure negative
mode does not force its selection as an ordinary AH step.

## A skipped negative root is not an unavailable direction

The split-only regression reached the high FeCl2 branch
($-2181.617645857189$ Eh) instead of the established lower branch. This is
not a successful performance result. It motivates checking the root filter,
not relaxing convergence conditions. The following rule uses only the
eigensystem that the regular AH iteration has already computed.

If the lowest augmented root is skipped by the existing reference-weight
filter, retain its unscaled tail $z$, images $H_0z,G_0z$, and reference weight
$|\alpha|$. Do not divide by a small or zero reference component. When
$z^TH_0z<0$, orient the ray so $g^Tz\leq0$ and let $t_{\max}$ place it on the
existing orbital/CI component boundary. Its quadratic change is

$$
\Delta m_z=t_{\max}g^Tz+\tfrac12t_{\max}^2z^TH_0z<0.
$$

Compare this value to the ordinary AH increment after its existing cap. Use
the ray only if it has strictly lower model energy. The normal start rule,
gradient-growth rollback, keyframes, physical trial energy and outer stopping
criteria remain. After taking a negative-curvature boundary step, the micro
loop ends and the accumulated trial enters the existing macro energy check.
Indeed, for $c=z^TH_0z<0$ and $g^Tz\leq0$, the updated slope is

$$
(g+tH_0z)^Tz=g^Tz+tc<0.
$$

Continuing to minimize the same indefinite quadratic would repeatedly prefer
the same ray; it has no finite unconstrained minimum. Ending that inner step
does not add a convergence test or a new physical evaluation beyond normal
macro acceptance.

No additional eigensolve, HVP, metric inverse, terminal
stability condition or tunable threshold is introduced.

The raw eigenpair residual belongs to $(\alpha,z)$, not to the boundary-scaled
step. A chosen boundary ray is explicitly marked nonconverged and reports
its complete unshifted model residual. Diagnostics record the skipped root,
curvature, both model changes and the direction choice. Those quantities,
plus the native FeCl2 trajectory, are needed to establish whether this is
actually the cause of the high-branch regression.

## Measured results and disposition

All new native runs used Hanhai25, 32 physical cores, the original singlet
input and the settings above. The existing versions were not rerun.

| Version | Job | Macrosteps | Orbital HVPs | Final energy / Eh | GNU wall / s | Peak RSS / KiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Historical horizontal-coordinate version, Hanhai22 | 1254797 | 17 | 242 | -2244.741071193477 | 1199.48 | 41439184 |
| Existing production baseline, Hanhai25 | 291577_0 | 21 | 316 | -2244.741075725727 | 759.22 | 41412664 |
| Positive-interior prototype | 291620 | 22 | 320 | -2244.741075726828 | 790.69 | 41413592 |
| Block-pure split images with cached negative-root ray | 291682 | 23 | 319 | -2244.741075727849 | 813.90 | 41412060 |

The 17-step row ran on EPYC 9654; Hanhai25 rows ran on EPYC 9755. Its wall
time is not a same-hardware performance reference. The 21-step baseline and
23-step trial both used `anode006`, but the original executable was built
with GCC 14.1 and the prototypes with GCC 14.3, with the same optimization
flags and the baseline Libcint/OpenBLAS libraries. Single-run times do not
isolate compiler or node-load effects. The step/HVP target was missed
independently of that timing limitation.

The final trial agrees with the baseline energy within $2.13\times10^{-9}$
Eh, but takes two more macrosteps, three more HVPs and 7.2% more wall time.
Through macrostep 18 it uses 224 HVPs versus the baseline's 255, but its energy
is still higher by $4.73\times10^{-6}$ Eh. Its remaining five macrosteps cost
95 HVPs versus 61 for the baseline's remaining three. Early action savings
therefore did not translate into faster convergence.

Of the final trial's 23 micro-sequence exits, 12 are keyframe mismatch, six
are micro budget, two are micro gradient, and one each is energy-model
mismatch, curvature boundary and gradient growth. The curvature ray is
selected once. These data do not isolate the separate causal effects of
splitting and the curvature policy; they show that the combined proposal
does not meet the requested performance target. Nor do they prove that a
tighter Ritz residual would fix the nonlinear model mismatch.

The final quick regressions (array 291678) all converged:

| Input | Macrosteps | HVPs | Final energy / Eh | GNU wall / s | Peak RSS / KiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| F2 | 3 | 19 | -198.751155830547 | 0.39 | 59380 |
| 241, repository guess | 3 | 20 | -230.720590533851 | 1.11 | 1017272 |
| 241, matched qualification guess | 4 | 45 | -230.720590533772 | 1.71 | 1039924 |
| FeCl2 | 14 | 192 | -2181.624464863246 | 3.79 | 440972 |

The earlier split-only regression ended on the high FeCl2 branch; the final
trial reached the established low branch without a new stability criterion.
That is a correctness qualification of this trial, not evidence of an iron
speedup. Seven algebra/transport/keyframe/HVP/diagonal tests passed before
the final native runs (build job 291672).

Raw outputs, action traces, hashes and resource logs are under
`results/diagnostics/ciah-subspace-20261017/`. The exact four changed source
files of the final trial are archived under `block-singlet/source/`, and
`block-singlet/source.sha256` identifies the full execution mirror. The final
trial binary SHA256 is
`b8a60e5a81fb7b59d2a218121c65f014bf61e51e3e2b8397541a0b119b490214`.
No failed experimental branch or compatibility switch remains in production.
Restoration was checked byte-for-byte against pre-experiment commit
`88af2d6e` for the four affected production/test files. The restored build
and the same seven unit tests passed on Hanhai25 (job 291696); no baseline
native performance run was repeated. The final trial executable is archived
alongside its logs as `block-singlet/xmvb-cpp.exe`.
