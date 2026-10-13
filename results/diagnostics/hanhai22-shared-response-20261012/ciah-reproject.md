# Reproject a changed CIAH RHS before expansion

Rejected experiment: branch `codex/ciah-reproject`, commit `6e764ad6`.
Based on 6bcc5bb4 (coordinate-norm CIAH), not a baseline rerun.
The main branch contains this report only; the experimental optimizer was not merged.
The nonorthogonal chart, HVP, generalized AH metric, preconditioner, clipping,
keyframe policy, root selection and external tolerances are unchanged.

## Single change

When a microstep or keyframe changes the gradient and an old residual is
pending, solve the projected AH problem using the new gradient and the retained
H/G images. Measure its complete residual, including outside the trial space.
Reuse the space if the original AH convergence criterion is satisfied, the
space is exhausted, or the residual is below the original dependence threshold.
Otherwise expand once using the refreshed residual and shift, not the old RHS.
No new tolerance or optional production mode is introduced.

Only one additional coordinate vector is retained per workspace, to identify
a changed RHS. No Hessian, metric, integral or per-thread large array is copied.
Projected solves may increase; expensive HVPs need not decrease. This is an
experiment, not a demonstrated performance improvement.

Tests use a nonidentity metric and an old residual outside the trial space.
One refreshed RHS is fully representable (zero new HVP); the other requires a
different outside direction (one HVP using the new residual). Dense generalized
AH roots, cached images and zero physical-norm callbacks are checked.

## Jobs and evidence

- Build and five unit-test executables: Hanhai22 job 1250595.
- Native parallel tests: array 1250597, index 4 FeCl2 and index 9 QOCGUM.
- Same inputs, 32 OpenMP cores, BLAS one thread, Libcint exact integrals,
  Davidson, analytic Hessian-diagonal preconditioner, gradient tolerance 1e-3,
  energy tolerance 1e-7 Eh. Source patch is in the frozen build directory.
- Reference: saved ciah-coordinate-1250538. No old version is run again.
- Unrelated local RI source changes are not included or transferred.

## Native result: reject for production

All five unit-test executables passed; both molecules converged to the correct
endpoint. Input and basis hashes match the saved reference runs. Timings are
single native runs on Hanhai22 AMD EPYC 9654 nodes, not repeated matched-node
measurements. RSS is GNU time's process maximum in MiB. No old version was run.

| Case | Steps, retained → experiment | HVP, retained → experiment | Wall / s, retained → experiment | RSS / MiB, retained → experiment | Experimental energy / Eh |
|---|---:|---:|---:|---:|---:|
| FeCl2 | 12 → 13 | 254 → 280 | 8.963 → 9.599 | 549.9 → 516.1 | -2181.624464869960 |
| QOCGUM | 4 → 4 | 52 → 51 | 15.230 → 17.346 | 1922.5 → 1971.2 | -480.094766681855 |

Energy differences are -2.322e-9 and +3.013e-12 Eh, respectively. The changed
expansion ordering preserves the tested endpoints but does not provide a useful
HVP reduction: FeCl2 requires an extra macroiteration and 26 more actions;
QOCGUM saves only one action. The timing measurements show no improvement,
although individual second-scale differences also include node/runtime noise.

Do not tune the reuse threshold to rescue this experiment. Keep 6bcc5bb4 as the
production algorithm; retain this implementation only on the experiment branch
for reproducibility, not as an optional backend or fallback.

Detailed numbers: [ciah-reproject-comparison.csv](ciah-reproject-comparison.csv).
Separate new-binary scalar audits: array 1250602, FeCl2 (4) and QOCGUM (15).
Their times are excluded from the native comparison.

## Reprojection audit

| New-binary diagnostic | RHS reprojections | Reused without expansion | Total projected solves | HVP |
|---|---:|---:|---:|---:|
| FeCl2 | 235 | 0 | 509 | 274 |
| QOCGUM | 39 | 0 | 90 | 51 |

The smallest reprojected residuals are about 9.57e-5 and 1.21e-4, respectively,
above the unchanged AH residual convergence threshold of 1e-6. The existing
space never passes the reuse test in these diagnostic runs. Additional small
solves therefore do not avoid operator actions. Refreshing expansion directions
changes the subsequent trajectory, but does not yield a useful measured gain.
FeCl2's diagnostic trajectory has 274 HVPs versus 280 in the independent native
run; audit counters are reported separately, not substituted into native data.

This rejects this conservative reprojection policy on these two cases, not all
forms of subspace reuse or a theorem that reuse cannot help. Production keeps
the validated 6bcc5bb4 optimizer. No benchmark of the old binary is submitted.
