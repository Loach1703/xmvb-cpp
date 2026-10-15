# FeCl2 AH tail: measured costs and global sample rank

## Latest bounded check: 1249390 / 1249395

The carrier-connection decomposition is now measured, not assumed; see
`FECL2_CARRIER_CONNECTION.md`. It reconstructs all 281 forcing vectors with
maximum coordinate error 2.14e-14. Basis re-representation is sizable but not
dominant, and the physical integral forcing remains comparable to the total.
No connection shear, dropped coupling or changed convergence rule was adopted.

Final-source scalar run 1249395 converged in 14 macros / 281 HVP, with energy
-2181.624464869471 Eh, GNU wall 10.55 s and RSS 512.3 MiB. Its complete structure
dual norms cost 0.8250 s (505 calls, 22,742 overlap actions). Small AH
spectrum/image combination cost 0.01712 s over 281 solves; projection/checks
cost 0.1492 s. The largest relative Galerkin residual was 2.6744e-12. AH
residual norm work, 0.6718 s, is included in the 0.8250 s all-norm total.

| Macro | Microsteps | Last predicted norm | Last AH residual norm | Orbital fraction of squared full residual |
|---|---:|---:|---:|---:|
| 11 | 20 | 1.77353e-3 | 1.77353e-3 | 96.68% |
| 12 | 21 | 2.53997e-3 | 2.53997e-3 | 99.50% |
| 13 | 25 | 4.92252e-4 | 4.92252e-4 | 98.50% |
| 14 | 26 | 1.34826e-4 | 1.34826e-4 | 97.71% |

These macros exit at the unchanged microstep budget, with a full residual still
above 9e-5. Accurate projected roots do not imply an accurate full-space Newton
step. The cost and residual evidence point to orbital complementary-space
expansion/preconditioning, not a need to tighten the tiny AH eigensolve or
replace the complete norm by a diagonal approximation. No speedup is claimed
for this diagnostic round; no old executable or iron-porphine long job was run.

## What the budget exit means

The installed PySCF `mcscf/newton_casscf.py::update_orb_ci` has a
gradient-dependent microstep budget and nonlinear gradient keyframes. Reaching
that budget commits the macro update; it does not certify a solved Newton
equation. `soscf/ciah.py::davidson_cc` retains trial vectors and Hessian images
and refreshes the gradient column of its projected AH problem. Our workspace
likewise retains its fixed-chart images. Orthogonal orbital rotation is not
copied into the nonorthogonal orbital update.

The native retained run 1249059 stops at keyframe mismatch in macros 9--11
and at the microstep budget in 12--13. The new run 1249116 has the same exit
sequence. Thus the tail is not exclusively a budget issue.

| Macro, 1249116 | Microsteps | AH basis | Source physical gradient norm | Last predicted norm | Last AH residual norm |
|---|---:|---:|---:|---:|---:|
| 12 | 23 | 25 | 1.90650e-3 | 5.79281e-4 | 5.79281e-4 |
| 13 | 26 | 28 | 5.83467e-4 | 2.42851e-4 | 2.42851e-4 |

Neither macro clips its steps. The last predicted norm remains above the
unchanged micro-gradient stop of 9e-5. The near equality to the AH residual
shows a substantial complementary-space residual, not an inaccurate small
eigensolve. The maximum relative Euclidean Galerkin residual in the whole run
is 4.40e-12. This comparison does not by itself prove that the complete
analytic Hessian is correct in every direction.

## Costs of the projected solve and physical inverse

Explicit wall-clock scopes in 1249116, not profiler sample percentages:

| Scope | Calls / actions | Wall time / s |
|---|---:|---:|
| AH projection and symmetry/rank checks | 257 | 0.1198 |
| Small generalized spectrum and image combination | 257 | 0.0140 |
| AH residual structure dual norm | 257 | 0.5395 |
| All structure dual norms | 625 | 0.9314 |
| Residual preconditioning | — | 0.1597 |
| Coupled HVP expansion, including sample ingestion | 257 | 2.0430 |

The residual-norm row is included in the all-norm row; do not add them.
The full structure norms performed 28,290 PCG iterations / 28,876 overlap
actions, with at most 53 iterations per solve. The largest certified true
relative residual was 1.48880e-8, below the unchanged square-root-of-machine-
epsilon target. These actions apply the overlap metric, not Hamiltonian sigma.
PySCF uses Euclidean norms in its orthonormal coordinates and therefore does
not have this inverse cost. That does not justify replacing the full metric
by the identity or the diagonal in our nonorthogonal chart.

The independently assembled bounded diagnostic in 1249025 previously found
no keyframe-decision disagreements with the matrix-free complete metric,
nor between the reported diagonal and complete metric decisions. Norm
inversion has measurable cost, but is not evidence explaining the excessive
HVP count. Fixing the small eigensolver alone cannot remove hundreds of HVPs.

## Failed inverse-scale experiment and the exposed rank defect

A positive RMS physical-length calibration of the base inverse introduced
no extra HVP. Runs 1249082, 1249089, 1249095 and 1249109 failed the existing
coarse projected-symmetry check. The local Gram-rank check did not resolve
the failure. In 1249109 the source/cache reciprocity discrepancy grew from
approximately 1e-11 to 1e-2, followed by a coarse skew of 27.30; the projected
roundoff bound was only 4.16e-10. It is not ordinary small-matrix roundoff.

The cache orthogonalizes directions and their previously evaluated images
with the same coefficients. A local independence check measures distance
from an already orthonormalized basis; it does not certify the conditioning
of the original sample collection from which the images were transformed.
Successive nearly dependent samples can amplify source-image errors through
the accumulated transformation.

Let normalized original sample directions be the columns of $D$, and let
$D=VR$ with orthonormal $V$. The complete sample Gram is

$$
D^T D=R^T R,\qquad
\frac{\lambda_{\min}(D^T D)}{\lambda_{\max}(D^T D)}
=\left[\frac{\sigma_{\min}(R)}{\sigma_{\max}(R)}\right]^2.
$$

The revised cache accumulates the small triangular factor $R$ during the
existing two-pass orthogonalization and checks this complete ratio against
the existing 1e-14 rank criterion. No HVP, physical direction, stopping
threshold or action precision is changed. Only redundant preconditioner
samples can be rejected; AH directions remain untouched. Capacity 60 adds
one shared 60-by-60 matrix (28.1 KiB), not a full-space dense metric or a
per-thread copy. An independent three-dimensional chained-sample test catches
a case that passes local independence but fails the complete Gram criterion.

With that check, the scale experiment 1249112 completed, but required 14
macros / 292 HVPs and was inferior to the retained optimizer. The scale class,
driver logic and unit test were deleted. No experiment switch remains.
The failed symmetry check was not widened or removed.

## Native comparison

All use original FeCl2 input, exact def2-TZVP integrals, Davidson, analytic
diagonal, 32 bound cores on Hanhai22, gradient tolerance 1e-3 and energy
tolerance 1e-7 Eh. Archived runs were not rerun. Different-node single timings
are exploratory; equal action counts are not an algorithmic speedup.

| Run | Macros | HVP | Energy / Eh | GNU wall / s | GNU RSS / MiB |
|---|---:|---:|---:|---:|---:|
| Retained 1249059, archived | 13 | 257 | -2181.624464862586 | 9.43 | 488.0 |
| RMS scale + global sample rank, 1249112, rejected | 14 | 292 | -2181.624464866758 | 10.50 | 486.9 |
| Original inverse + global sample rank, 1249116 | 13 | 257 | -2181.624464862153 | 9.45 | 541.3 |

1249116 passed five compute-node unit executables, including the global
sample-rank regression. Its energy differs from the independent PySCF
coupled-Newton endpoint -2181.624464883733 Eh by +2.158e-8 Eh.
PySCF default-AH endpoint -2181.624464702896 Eh is also retained as a separate
reference; these two PySCF endpoints themselves differ by 1.808e-7 Eh.

The frozen 1249116 executable also passed the parallel small-system
regression array 1249123 on Hanhai22:

| Case | Macros | HVP | Energy / Eh | GNU wall / s | GNU RSS / MiB |
|---|---:|---:|---:|---:|---:|
| F2 | 3 | 23 | -198.751155830526 | 0.28 | 17.2 |
| 241 | 3 | 21 | -230.720590392836 | 2.65 | 809.2 |
| MnF2 | 10 | 137 | -1348.893353223951 | 6.59 | 505.5 |

All converge with exit zero and retain the macro/HVP counts and energies of
the archived 1249065 regressions. Times/RSS are separate single-node native
observations, not proof of a speedup or a statistically established memory
change.

## Conclusion and remaining work

### Rejected complementary-space inverse experiment

Follow-up 1249133 tested the exact constrained inverse of the same
approximate shifted operator. For orthonormal trial basis $Q$, inverse
preconditioner $P$ and residual $r$, the approximate KKT equations are

$$
P^{-1}t+Q\mu=r,\qquad Q^Tt=0,
$$

giving

$$
t=Pr-PQ(Q^TPQ)^{-1}Q^TPr.
$$

This is generally different from orthogonalizing $Pr$ afterwards. It is a
search-subspace-restricted correction, not a complete inner Jacobi--Davidson
solve. The motivation is discussed in Genseberger and Sleijpen,
[Alternative correction equations in the Jacobi-Davidson method](https://ir.cwi.nl/pub/4645).
The physical AH operator, nonorthogonal metric, retraction, keyframes and
all stopping rules remained unchanged. Only the expansion preconditioner
changed. Each expansion still made one HVP, but applied the existing inverse
to its RHS and to every cached basis vector.

An independent dense augmented KKT solve verified the new direction. The
five compute-node unit executables passed. The first batch 1249131 stopped
at a legacy unit assertion expecting only one preconditioner call; it did
not run chemistry. The corrected accounting assertion in 1249133 requires
one HVP and the mathematically necessary inverse calls.

| Run | Macros | HVP | Energy / Eh | GNU wall / s | GNU RSS / MiB |
|---|---:|---:|---:|---:|---:|
| Archived 1249116 | 13 | 257 | -2181.624464862153 | 9.45 | 541.3 |
| Restricted inverse 1249133, rejected | 25 | 362 | -2181.624464701934 | 16.03 | 500.7 |

Preconditioning cost rose from 0.160 to 1.817 s; even removing that added
cost would not justify the increased macro/HVP counts. The experiment is
not retained: its implementation and unit-test modifications were deleted
locally and from the execution mirror, with no switch or fallback. The
default inverse is unchanged. This result rejects this specific constrained
expansion as a FeCl2 improvement; it does not rule out all Jacobi--Davidson
correction solvers or establish a derivative defect.

### Recycled structure direction: missing root phase

The unchanged-model diagnostic 1249162 found that the relaxed CI root was
opposite in phase to the last coupled coefficients in macros 2--6. The
complete physical overlap with the reference was approximately -0.98 to
-1.00, not an ambiguous near-zero overlap. The driver recycled the previous
raw CI direction without reversing it, while recycling the orbital direction
unchanged. This does not affect the objective or derivatives, but makes the
mixed initial direction depend on the eigensolver's arbitrary eigenvector
phase.

Before projecting the next CI direction into the new tangent chart, it must
transform with the relaxed root. In a fixed raw coefficient representation,

$$
C_{\mathrm{new}}=\sigma\widetilde C,\qquad
\delta C_{\mathrm{new}}=\sigma\delta\widetilde C,\qquad
B_{\mathrm{new}}=\sigma\widetilde B,
$$

where $\sigma$ is -1 for negative reference/root overlap and +1 otherwise.
Consequently $B_{\mathrm{new}}^T\delta C_{\mathrm{new}}$ is phase invariant.
The driver now orients the recycled response using
$C_{\mathrm{ref}}^T M_{\mathrm{new}}C_{\mathrm{new}}$. It reads the existing
Davidson overlap/root image; dense reference solves use their already owned
overlap matrix. No HVP, integral transformation, dense assembly or metric
inverse is added. The reference uses the actual accepted displacement,
which already contains any line-search scaling.

Unit tests verify positive and negative root phases, an unchanged orbital
direction and covariance of the coupled orbital image. Five unit executables
passed in 1249185. The matrix-free FeCl2 result was 14 macros / 282 HVPs,
E=-2181.624464870537 Eh, GNU wall 9.99 s, RSS 520.1 MiB. That is not a speed
improvement over 13 / 257. The phase fix is a restart-consistency correction,
not a claimed explanation or solution of slow tail convergence.

The new dense reference path 1249211 also converged: 12 macros / 244 HVPs,
E=-2181.624464864200 Eh, GNU wall 13.49 s, RSS 533.0 MiB. Its difference from
the new Davidson result is 6.337e-9 Eh. This is a new-code path regression,
not a rerun of an old performance baseline.

1249208 tested the phase-corrected binary on F2, 241 and MnF2. All converged:
3 / 23, 3 / 21 and 10 / 144 macro/HVP counts, respectively. Energies remain
-198.751155830527, -230.720590392836 and -1348.893353223965 Eh. The MnF2
HVP count increased from the archived 137, so no universal speedup is claimed.
This binary's phase reference was subsequently corrected to avoid applying
the line-search scale twice; 1249221 verifies the final implementation.

Final-source 1249221 passed all five unit executables and converged FeCl2
in 14 macros / 281 HVPs, E=-2181.624464869146 Eh, GNU wall 9.93 s, RSS
510.0 MiB. Its difference from the retained independent PySCF coupled-Newton
reference is +1.4587e-8 Eh. The final-source scaled-line-search regression
1249224_2 converged MnF2 in 10 macros / 144 HVPs, E=-1348.893353223967 Eh,
GNU wall 5.84 s, RSS 498.1 MiB. These are correctness/trajectory observations,
not a claimed speedup; archived baselines were not rerun.

An additional controlled ablation 1249203 omitted only the recycled CI
increment, while retaining the complete coupled operator and the recycled
orbital direction. It required 16 macros / 306 HVPs, E=-2181.624464800146 Eh,
GNU wall 11.08 s, RSS 499.1 MiB. It was rejected and removed; the phase-corrected
CI direction remains. No orbital-only bootstrap switch was introduced.

This round fixes certification of cached-image sample conditioning, not the
slow tail contraction. The micro-budget exit remains; increasing budgets,
relaxing norms, or accepting the ineffective RMS scale is not a solution.
The evidence prioritizes complementary orbital/coupled-space preconditioning
and frozen-curvature/keyframe consistency, rather than LAPACK tuning or
metric-tolerance changes. Earlier failed RHS-refresh and split-AH trials
must also remain part of that assessment, not be repeated as new ideas.

No new iron-porphine long run was launched: the required FeCl2 improvement
has not been demonstrated. Raw logs, hashes, build flags and GNU times are
retained in the numbered result directories. `summarize_ritz_cost.py`
reproduces the scalar accounting without manually transcribing iterations.

## Auxiliary shift experiment and exact norm reuse

1249250 tested replacing the auxiliary identity shift in the full-OEO
analytic-diagonal inverse by the physical flag-metric shift. In orthonormal
frame components, core--active motion has physical weight two, while
core--external and active--external motion have weight one. A bounded dense
unit test verified the shifted model inverse against the independently
applied physical metric. All five unit executables passed.

This is a change to an approximate preconditioner, not a correction to the
complete AH equation: Hessian actions, the generalized Ritz metric and
physical residual norms already use the declared physical metric. An
auxiliary identity-shifted preconditioner is a legitimate approximation.
The more literal shift was inferior and was removed, together with its
experimental tests; no switch was retained.

| Run | Macros | HVP | Energy / Eh | GNU wall / s | GNU RSS / MiB |
|---|---:|---:|---:|---:|---:|
| Archived phase-correct 1249221 | 14 | 281 | -2181.624464869146 | 9.93 | 510.0 |
| Physical flag shift 1249250, rejected | 17 | 305 | -2181.624464865874 | 11.47 | 489.9 |

The new performance change uses an exact identity rather than replacing a
norm or loosening an inverse tolerance. For a generalized AH candidate,

$$
r=\alpha(g+Hp-\omega Gp),\qquad
g_+=g+Hp=\frac{r}{\alpha}+\omega Gp.
$$

Here $\alpha$ is the signed reference coefficient, $\omega$ the AH Ritz
eigenvalue, $G$ the complete coupled physical metric, and $p$ the unscaled
candidate. Therefore

$$
\|g_+\|_{G^{-1}}^2=
\frac{\|r\|_{G^{-1}}^2}{\alpha^2}
+\frac{2\omega}{\alpha}r^Tp
+\omega^2 p^TGp.
$$

The residual norm and $Gp$ are already available. The cross term is retained
instead of assuming exact Galerkin orthogonality. This removes one repeated
complete metric inverse per unscaled accepted microstep, with scalar/vector
arithmetic only and no extra matrix storage. A clipped candidate is a
different affine gradient and continues to require its own complete norm.
The independently solved dense metric in the CIAH unit tests checks the
identity for anisotropic metrics, indefinite curvature and refreshed RHSs.

This optimization does not claim to remove the complementary-space residual
or reduce the HVP count. PySCF retains its gradient-dependent microstep
budget and cached Davidson images; see
[update_orb_ci](https://pyscf.org/_modules/pyscf/mcscf/newton_casscf.html)
and [davidson_cc](https://pyscf.org/_modules/pyscf/soscf/ciah.html).
The extra inverse is specific to our nonorthogonal metric, not a missing
PySCF eigensolver tolerance.

### Compute-node validation of norm reuse

1249298 passed all five unit executables. The original-guess FeCl2 native
run converged in 14 macros / 282 HVPs, E=-2181.624464869792 Eh, GNU wall
9.97 s, RSS 491.6 MiB. Relative to archived 1249221, complete structure
norm calls decreased from 682 to 506, overlap actions from 31,418 to
22,791, and scoped inverse time from 1.0167 to 0.7531 s. The maximum true
relative solve residual remains 1.48971e-8; no precision was relaxed.
Total time is not lower, so this is demonstrated removal of redundant
inverse work, not demonstrated overall speedup or reduced HVP count.

The frozen executable regression 1249320 completed all three cases:

| Case | Macros | HVP | Energy / Eh | GNU wall / s | GNU RSS / MiB |
|---|---:|---:|---:|---:|---:|
| F2 | 3 | 23 | -198.751155830526 | 0.29 | 17.0 |
| 241 | 3 | 21 | -230.720590392835 | 2.63 | 815.0 |
| MnF2 | 10 | 144 | -1348.893353223966 | 5.50 | 547.2 |

Counts and endpoints agree with the retained phase-correct regression.
No archived executable was rerun.

1249325 adds physical orbital/structure norm accounting without an extra
inverse or HVP. Five units again passed; FeCl2 converged in 13 / 254,
E=-2181.624464846043 Eh, GNU wall 9.10 s, RSS 514.8 MiB. Its difference
from the independent coupled PySCF endpoint is +3.7690e-8 Eh. Do not
attribute its lower HVP count to logging or to norm reuse: the first norm-
reuse run had 14 / 282, and no new direction-construction policy was added.
The diagnostic explicitly sequences the two norm evaluations and records
their scalar components; this does not change their mathematical definition.

The final AH residual blocks in the instrumented run are:

| Macro | Orbital dual norm | Structure dual norm | Orbital fraction of residual squared |
|---|---:|---:|---:|
| 11 | 1.96029e-3 | 2.21482e-4 | 98.74% |
| 12 | 2.54112e-3 | 1.91620e-4 | 99.43% |
| 13 | 4.01252e-4 | 4.29475e-5 | 98.87% |

Its maximum relative Galerkin residual is 3.0581e-12. The last two macros
still stop at the microstep budget, with last predicted dual norms
approximately equal to the complete AH residual norms. This localizes the
unresolved residual primarily to the **orbital complement of the AH trial
space**, not an unsolved structure equation or inaccurate Ritz eigenpair.
It does not establish whether a particular improved orbital preconditioner
or nonlinear model update will cure it.

Scoped costs in 1249325 are 0.1185 s for AH projection/checks, 0.0137 s
for spectrum/image combination, 0.6903 s for all structure norm inverses,
and 1.9967 s for coupled expansions including sample ingestion. Chart
construction is separately 1.8983 s; its 0.4240 s positive-block factor
scope is nested, not additive. These are explicit wall-clock scopes,
not profiler percentages. Removing unused chart factors would be a cost
optimization only and requires checking the independent block-LBFGS APIs;
it cannot explain or remove the complementary residual.

No iron-porphine long run is started on the basis of these modest/variable
FeCl2 timings. The requested robust reduction in HVPs remains outstanding.

## Orbital-block inverse experiment

Two short, original-guess tests replaced only the analytic orbital
preconditioner with a signed orbital-block inverse; the actual coupled AH
operator, metric, retraction, keyframes and thresholds stayed unchanged.
Local spectral factors already computed by the chart were retained, rather
than recomputing a Hessian or creating a dense global orbital matrix.

1249343 used the existing local Euclidean shift. This changes both the
block approximation and the auxiliary shift representation, and therefore
is not an isolated test of block Jacobi. 1249346 instead tested genuine
block Jacobi of the same full-frame shifted model. With frame map $T$,
diagonal frame curvature $D$ and per-orbital column block $T_i$, its inverse
uses

$$
K_i=T_i^TDT_i,\qquad J_i=T_i^TT_i,\qquad
P_i(\omega)=(K_i-\omega J_i)^{-1}.
$$

A generalized local eigendecomposition satisfies

$$
K_iW_i=J_iW_i\Lambda_i,\quad W_i^TJ_iW_i=I,\quad
P_i(\omega)=W_i(\Lambda_i-\omega I)^{-1}W_i^T.
$$

The independent bounded dense pullback test was updated to check exactly
this inverse, while retaining the original chart and component-cap tests.
Five unit executables passed in both runs.

| Run | Macros | HVP | Energy / Eh | GNU wall / s | GNU RSS / MiB |
|---|---:|---:|---:|---:|---:|
| 1249343, local Euclidean shift | 15 | 325 | -2181.624464010912 | 10.63 | 499.6 |
| 1249346, same-model block Jacobi | 16 | 346 | -2181.624464648738 | 12.53 | 489.2 |

Neither improves on retained 13--14 macro / 254--282 HVP observations.
The first endpoint differs from independent PySCF coupled Newton by
+8.72821e-7 Eh; the second by +2.34995e-7 Eh. Against the separately
retained PySCF default-AH endpoint, the differences are +6.91984e-7 and
+5.41577e-8 Eh, respectively. These are optimizer endpoint differences,
not a change in Hamiltonian or spin space. A completed run under the shared
criteria is not proof of a particular absolute-energy error.

Both experimental implementations and their tests were removed from the
workspace and execution mirror. No block-Jacobi selector, parameter or
fallback is retained. The result supports keeping the current full-frame
inverse rather than dropping its cross-orbital representation merely for
local simplicity; it does not rule out other better orbital models.

## Metric-orthogonal AH experiment

1249350 replaced Euclidean basis orthogonalization by two-pass physical
$G$ orthogonalization, updating cached metric images with the same
coefficients. Five units passed, including independently checked physical
normalization/orthogonality. The actual FeCl2 trial stopped before accepting
a macro, after eleven HVP expansions, at the existing projected-metric
symmetry check. This exposes cancellation amplification in the experimental
image recurrence; it is not evidence that the retained metric operator is
asymmetric. The check was not widened.

1249359 reapplied the actual metric to each admitted vector before its HVP,
rather than trusting its recursively subtracted image. This adds a second
metric action per admitted direction, but no extra HVP. Five units passed
and FeCl2 converged in 14 macros / 281 HVPs, E=-2181.624464863864 Eh,
GNU wall 9.92 s, RSS 496.5 MiB. The final four macros still stop at the
microstep budget. It provides no HVP reduction over the retained algorithm;
its extra metric work is therefore not retained. Both orthogonalization
experiments and their test changes were removed locally and from the mirror.

These experiments narrow the next action: compare physical integral coupling
with representation-induced carrier coupling before changing another
preconditioner or eigensolver policy. The exact decomposition and selected-
space caveats are recorded in `FECL2_CARRIER_CONNECTION.md`; they are not
yet a validated diagnosis or an implementation change.
