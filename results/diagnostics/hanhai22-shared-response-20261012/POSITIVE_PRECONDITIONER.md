# Positive one-electron CIAH preconditioner restoration

## Scope

Production commit `ab7e3c09` changes only the explicit `one-electron`
orbital residual preconditioner in the single-state CIAH backend. It restores
`OrbitalChart::apply_inverse_reduced_block_preconditioner`. The analytic
Hessian-diagonal default, CI preconditioner, full HVP, chart, orbital update,
CIAH policies and convergence tolerances are unchanged. No molecular rule or
new optimizer parameter is introduced. No large array is copied per thread.

The existing positive block is

$$
P_i=V_i\operatorname{diag}\!\left(\max(|\kappa_{ij}|,\tau_i)\right)V_i^T,
\qquad
\tau_i=\sqrt{\epsilon_{\mathrm{mach}}}\max_j|\kappa_{ij}|.
$$

Here $\kappa_{ij}$ and $V_i$ are the local one-electron model eigenvalues and
eigenvectors. The restored orbital callback returns $P^{-1}r$, rather than
inverting the signed local spectrum shifted by the AH eigenvalue. This does
**not** make the actual Hessian positive.

## Separate metric checks

The complete coupled AH problem and residual still contain the physical metric:

$$
\begin{pmatrix}0&g^T\\g&H\end{pmatrix}
\begin{pmatrix}\alpha\\v\end{pmatrix}
=\omega
\begin{pmatrix}1&0\\0&G\end{pmatrix}
\begin{pmatrix}\alpha\\v\end{pmatrix},
\qquad r=\alpha(g+Hp-\omega Gp),\quad p=v/\alpha.
$$

The displacement norm uses $p^TGp$; the residual norm uses the physical dual
metric. Neither is replaced by the preconditioner or by the identity matrix.
The existing `test_ciah` independently checks the generalized root, metric
image and residual against dense references with an anisotropic metric.

For a shifted positive preconditioner, the relevant model is

$$
(P+\sigma G)d=r.
$$

`shifted_block_inverse` alone instead inverts $P+\sigma I$. It is suitable
as the PCG block preconditioner, not as an exact inverse of the full shifted
physical-metric model. The added `check_physical_metric_shift` test uses an
off-diagonal, non-unit metric: it matches a dense solve of $P+\sigma G$
and explicitly excludes the $P+\sigma I$ result. Existing rescaling and
zero-shift tests are retained.

An isolated diagnostic patch, `ciah-metric.patch`, tests
$\sigma=\max(0,10^{-8}-\omega)$ with the complete metric. The $10^{-8}$ is
the existing CIAH level shift, not a new policy. This experiment requests
machine epsilon as the inner budget; the existing solver consequently uses
its machine-derived square-root-epsilon residual floor. It is an attribution
experiment, not a production accuracy change or recommended inner policy.
The patch is applied only to the remote execution mirror for compilation and
automatically reversed. Production remains the unshifted positive interface.

## Execution and validation

- Cluster: Hanhai22; Slurm partition `CPU-192C768GB`.
- Build/tests: `1247629`; CIAH, sparse orbital quotient and shifted-metric
  preconditioner tests all passed.
- Positive block regressions: `1247631`, five systems in parallel.
- Complete-metric diagnostic build: `1247645`; runs: `1247646` (MnF2, FeCl2).
- Earlier metric build `1247642` failed during patch dry-run, before compilation
  or chemistry; its blocked dependent job was cancelled and resubmitted.
- Same frozen inputs as `1247579`, exact integrals, 32 allocated physical cores,
  explicit core binding, OpenMP 32 and BLAS 1. Outer thresholds remain
  $10^{-3}$ and $10^{-7}\ E_h$.
- Historical baseline `1247579` is reused, not rerun. Nodes differ, so single-run
  wall times are not a strict same-node speed ratio. HVP counts are unaffected
  by that qualification.
- Binary/input hashes, build flags, nodes, thread settings, native GNU-time
  output and raw NEO traces are retained in the corresponding artifact folders.
- Unrelated local RI changes were not transferred or committed.

## Restored positive-block results

| System | Steps | Orbital HVP | Energy (Eh) | Wall (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|
| F2 | 3 | 23 | -198.751155830526 | 0.28 | 17.17 |
| 241 | 3 | 23 | -230.720590392889 | 2.31 | 818.66 |
| 240 | 4 | 37 | -343.446302353107 | 6.53 | 1530.35 |
| MnF2 | 18 | 288 | -1348.893353118897 | 9.19 | 478.85 |
| FeCl2 | 23 | 574 | -2181.617644476738 | 16.73 | 423.37 |

All satisfy the existing stopping criteria. **FeCl2 does not reach the correct
low-energy branch**, despite the nominal convergence report: it is about
6.820 mEh above the existing analytic-diagonal CIAH endpoint. This is not a
successful FeCl2 energy correction. Restoring positivity resolves the extreme
iteration stagnation but is insufficient to guarantee the same stationary point.

Existing analytic-diagonal results: MnF2 11 steps/169 HVP, 5.56 s,
551.84 MiB, -1348.893353223710 Eh; FeCl2 20 steps/429 HVP, 11.91 s,
385.86 MiB, -2181.624464318884 Eh. They remain the correct reference here.
No extra stability/convergence condition has been added.

## Complete-metric shift experiment

| System | Steps | Orbital HVP | Energy (Eh) | Wall (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|
| MnF2 | 21 | 351 | -1348.893353221788 | 11.91 | 459.37 |
| FeCl2 | 17 | 418 | -2181.617644026137 | 13.81 | 433.91 |

The full-metric shifted positive model does not restore the low-energy FeCl2
branch either. Its lower FeCl2 action count cannot be called an improvement
over the analytic-diagonal reference, because the endpoints differ by about
6.820 mEh. It is also not faster than the unshifted positive model on MnF2.

This separates a confirmed signed/positive interface regression from an
unresolved root/direction-coverage problem. The two positive experiments do
not establish that the metric or the exact HVP is incorrect. They also do not
prove that a positive preconditioner must always find a particular minimum.
The exact HVP's physical negative curvature is preserved in both experiments.
Further diagnosis should examine the actual orbital directions entering the
AH subspace and their overlap with the directions that lead to the lower
FeCl2 solution; it must not claim an unmeasured signed spectral pole or add
extra convergence conditions to conceal the differing stationary points.

The unshifted/shifted runs jointly change whether a shift is applied, so their
performance difference alone is not an isolated test of $G$ versus $I$.
That distinction is established only by the independent dense-reference unit
test in this round. A shifted-identity molecular ablation was not run.
