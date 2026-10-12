# FeCl2 AH tail: measured costs and global sample rank

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
