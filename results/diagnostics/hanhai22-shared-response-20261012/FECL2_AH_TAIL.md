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
