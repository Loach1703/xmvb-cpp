# One-electron preconditioner regression: code and model audit

## Confirmed change, not a parameter diagnosis

The CIAH integration in `a8d60c8c` changed the meaning of `one-electron`.
The current paired run is not a test of the historical positive one-electron
block preconditioner inside the new CIAH workflow. Production optimizer code
and tolerances are not modified by this audit.

Historical references inspected: `c2776c4` (saved paired CSV) and `5568b9c6`
(latest prior correct-energy NEO). Both construct a positive orbital model
from the absolute spectrum of local one-electron curvature. They use the
block inverse at zero shift. On the coupled response path, the shifted model
also includes the full physical orbital metric through a matrix-free inner
solve; the local identity-shifted block inverse only preconditions that solve.

## Mathematical distinction

Let $K_i=V_i\operatorname{diag}(\kappa_{ij})V_i^T$ be the symmetric local
one-electron surrogate curvature pulled back to the orbital chart. It is not
the complete orbital Hessian: the two-electron and orbital/structure response
contributions need not have the same spectrum.

The historical positive block is

$$
P_i=V_i\operatorname{diag}\left[\max(|\kappa_{ij}|,\tau_i)\right]V_i^T,
\qquad
\tau_i=\sqrt{\epsilon_{\mathrm{mach}}}\max_j|\kappa_{ij}|.
$$

With $P=\operatorname{blockdiag}(P_i)$, orbital residual $r$ and physical
orbital metric $G$, the historical unshifted direction is $P^{-1}r$. The
coupled solver's shifted preconditioning approximately solves

$$
(P+\sigma G)d=r,\qquad \sigma\geq0.
$$

The present CIAH call instead uses, independently for each local block,

$$
d_i=V_i\operatorname{diag}\left[
f(\kappa_{ij}-\omega+\delta)^{-1}
\right]V_i^T r_i,
\qquad \delta=10^{-8},
$$

where $\omega$ is the projected augmented eigenvalue and $f(t)=t$ except
that entries with $|t|<10^{-8}$ are replaced by $10^{-8}$. The code retained
the positive block data but overwrote the AH signed spectral fields with the
one-electron surrogate's signed eigenvalues, then selected the AH inverse for
every non-identity preconditioner.

Thus this is an actual signed surrogate inverse, not the former positive
one-electron inverse. The two can differ substantially when the surrogate
has negative modes or eigenvalues close to the AH shift. Unlike $P$, the
signed shifted surrogate need not be positive or safely separated from zero.
The current projected AH problem still uses the full physical metric: there
is no evidence from this code audit that the exact Hessian, orbital chart,
gradient or objective was changed by choosing the preconditioner. However,
the shifted *preconditioner* no longer includes global $G$ coupling. These
are separate statements.

For a true lowest generalized AH eigenvalue, the corresponding Hessian shift
has useful positivity properties. Those properties do not automatically hold
for an unrelated one-electron surrogate, nor for a higher root selected by
the reference-weight filter. A signed Davidson preconditioner is not
mathematically forbidden; the confirmed defect is the unvalidated semantic
replacement of the proven positive one-electron policy. This audit does not
yet prove near-pole amplification in either molecule: that needs actual
spectral/action measurements or controlled ablation.

## Trace evidence

Same frozen binary/input/32-core allocation, paired job 1247579. Only the
preconditioner option differs within each pair. Initial physical gradient
norms agree between variants to the recorded floating-point precision.

| FeCl2 first macro | Physical step norm | Actual decrease (Eh) | Actual/predicted decrease | End physical gradient |
|---|---:|---:|---:|---:|
| Analytic Hessian-diagonal | 0.3700808662 | 2.4059289444e-5 | 0.99238109 | 0.0049903228 |
| Current signed one-electron | 0.0008460543 | 2.2822127903e-6 | 1.00500424 | 0.0053377386 |

The step norms differ by about 437x. Both local energy predictions are
accurate; the signed one-electron run has no rejected macro trials through
the first 725 recorded steps, yet little useful progress. The data do not
support blaming excessive energy rejection or inaccurate forward energy.
They support investigating subspace directions and coverage instead.

MnF2 also has slow progress rather than repeated rejection: through 800
recorded steps it has only three rejected trials, with physical gradient
still around 0.017. The analytic variant converges in 11 steps with final
physical gradient about 8.6e-5. These are running snapshots, not converged
one-electron results.

## Which old results are usable?

The saved `results/history/paired_runs.csv` records old NEO MnF2 in 11 steps
and 199 orbital HVPs, at -1348.893351161935 Eh. This is a useful algorithmic
reference, although current source/cluster/trajectory are not identical.

Its old FeCl2 result (4 steps, 72 HVPs, -2181.617639896249 Eh) is on the
wrong high-energy branch, approximately 6.8 mEh above the reference. It
cannot be used as a correct-energy performance target. The prior repaired
NEO job 1247396 reaches -2181.624464422983 Eh in 17 steps, but at high cost.

## Clean next verification

Keep the new CIAH action, root selection, chart, micro/keyframe policy, CI
preconditioner, initial guess and outer tolerances fixed. First compare the
current signed one-electron inverse with the existing historical positive
block inverse on identical residuals at the same accepted point. Record
residual/vector pairing, physical direction norm and shifted-model residual.
Then change only the orbital preconditioning callback in a diagnostic
ablation, testing positive block and its full-metric shifted counterpart.

This separates the sign change from the omission of global metric coupling
and from the CIAH framework change. It is not a molecule-specific policy,
threshold search, old-version rerun or a claim that restoration must solve
every convergence problem. No production repair is made before this test.
