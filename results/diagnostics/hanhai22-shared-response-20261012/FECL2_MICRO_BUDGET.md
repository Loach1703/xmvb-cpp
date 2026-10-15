# FeCl2: micro-budget exits versus AH subspace accuracy

## Scope

Use the archived native runs 1249512 and 1249543. Do not rerun the old
implementation. The new timing-only short run is 1249566 on Hanhai22; it
rebuilds and tests the current source before running from the original guess.
No numerical policy, tolerance, HVP, coordinate map or physical update changes.
The additional events allocate no arrays and evaluate no additional operators.

Reference code: PySCF 2.14.0 `mcscf/newton_casscf.py:update_orb_ci` and
`soscf/ciah.py:davidson_cc`. Sources were read from the installed reference
environment on Hanhai22, including the generator's post-yield control flow.

## The exit is not a converged Newton solve

In 1249543 macro 11, there are 24 applied microsteps and 26 admitted AH
directions. The maximum AH dimension remains 30. After the second keyframe,
the keyframe gradient norm is 0.0006044028549104222. The unchanged PySCF
micro-budget formula gives

$$
m_{\max}=\max\{10,10-\operatorname{int}[2\log(\|g_{\mathrm{KF}}\|+10^{-7})]\}=24.
$$

The final predicted full norm is 0.0003134139734294, still above the
unchanged micro-gradient target 0.00009. Thus the budget, not the AH
dimension limit or the micro-gradient criterion, ends this macro.
This is an approximate CIAH macrostep, not a fully solved Newton equation.
Successful outer convergence does not certify an exact inner Newton step.
Conversely, a budget exit alone is not proof of incorrect derivatives.

## Projected and full residuals are different

For physical metric $G$, current gradient $g$, Ritz step $p$ and AH root
$\omega$, the residual is

$$
r=\alpha(g+Hp-\omega Gp).
$$

The small generalized eigensolve enforces $Q^Tr\simeq0$. This does not imply
$\|r\|_{G^{-1}}\simeq0$ outside the current trial space. Across 1249543 the
maximum relative Euclidean Galerkin residual is 4.11e-12. The final full
residual is 0.0003134139661297, dominated by its orbital part:

| Component | Complete physical dual norm |
|---|---:|
| Orbital | 0.0003127431791414 |
| Structure | 0.0000204943422852 |

Orbital contribution is 99.5724% of the squared residual norm. The final
residual is about 40.7% of the macro's initial full gradient norm, rather
than a rapidly vanishing inexact-Newton forcing ratio.

The second macro-11 keyframe changes the gradient by 0.0004970544503074.
Only 0.0001042850453852 is captured in the existing physical AH trial space:
95.5981% of this correction's squared norm lies outside it. This exposes
both slow orbital-space enrichment and a changing right-hand side. It is
not evidence that the small projected eigenproblem was solved inaccurately.
It also does not prove the immutable origin Hessian is erroneous: a real
nonlinear keyframe legitimately changes the gradient and local curvature.

## Native cost accounting

Timers overlap where explicitly noted; do not sum nested measurements.

| Measured phase / s | Archived raw CI inverse, 1249512 | Carrier CI inverse, 1249543 |
|---|---:|---:|
| Complete structure dual norms | 0.1131 | 0.1228 |
| All Ritz norm sections, including those norms | 0.2591 | 0.2372 |
| Incremental tall-vector projection updates | 0.0317 | 0.0286 |
| Ritz projection checks | 0.0092 | 0.0076 |
| Small spectrum plus Ritz-image construction | 0.0234 | 0.0185 |
| AH inverse applications | 0.2173 | 0.2917 |
| Expansion HVPs, including image observer | 2.5089 | 2.0858 |
| Accepted-point analytic diagonal builds | 1.9359 | 1.7033 |
| Native whole process | 9.92 | 8.76 |

All structure norm solves need at most one iteration; the largest certified
relative residual in 1249543 is 9.79e-16. No norm precision is relaxed.
Eliminating every Ritz norm section could save at most about 2.7% of this
native run. Eliminating the small spectrum/image section alone could save
at most about 0.21%. These are measured-section upper bounds, not predicted
achievable speedups or exclusive hotspot shares.

## PySCF comparison and next evidence

The current budget formula, budget-before-keyframe ordering, retained
`xs`/`ax` equivalent and delayed residual expansion match PySCF's coupled
Newton implementation. PySCF also permits approximate budget-limited
macrosteps. Our physical metric and nonorthogonal chart mean that equal
policy constants do not imply equal trajectories or residual components.

Do not increase the budget, weaken the norm or replace the full norm with
a diagonal one. The remaining mathematical target is better enrichment of
the orbital residual outside the current space, including after keyframe
RHS changes. Previous rejected RHS-refresh/rank-two experiments must not
be relabeled as new solutions.

1249566 completed on anode31: all five unit executables passed; application
exit zero, 13 macros/252 HVPs, energy -2181.624464872092 Eh, GNU wall 9.13 s,
peak RSS 504.66 MiB. It logs the actual budget, basis dimension, micro target
and AH convergence flag at each exit. The final macro exits at `micro_gradient`
(8.3474e-5 < 9e-5), not the budget; macro 12 exits at budget 25 with residual
2.6995e-4. Two earlier macros end through gradient growth instead.

Across all complete dual-norm calls, orbital work costs 0.2721 s and structure
work 0.1010 s. These cover more than the Ritz-only norm timer (0.2225 s), which
must not be added again. Spectrum/image work costs 0.0216 s.

The timer-only rebuild produces a different trajectory from 1249543, despite
unchanged numerical statements and policies. Its source diff contains clocks
and scalar output only. The source-level invariance is not proof of bitwise
parallel reproducibility; no attribution to rounding/scheduling has been
validated. Do not count the 11-to-13 macro difference as a mathematical change
or claim a universally repeatable eleven-macro result. The endpoint remains
1.1641e-8 Eh above the PySCF coupled-Newton reference.

Existing iron jobs 1249554 remain bounded four-step diagnostics, not convergence
claims. Singlet completed four macros/51 HVPs at -2244.719899599181 Eh, 433.25 s,
38.6219 GiB GNU RSS; triplet was still running at the last scheduler check.
The archived singlet four-step run has 36 HVPs/424.92 s on a different node.
Lower fourth-step energy is not a whole-run speedup or convergence proof.
