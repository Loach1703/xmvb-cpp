# CIAH coordinate-norm control

## Scope and mathematical basis

This change removes physical dual-norm inversion from single-state CIAH
microiteration control. It does not replace the nonorthogonal metric with the
identity in the augmented Hessian equation. In the origin chart, that equation is

\[
\begin{pmatrix}0&g^T\\g&H\end{pmatrix}
\begin{pmatrix}\alpha\\\alpha p\end{pmatrix}
=\omega\begin{pmatrix}1&0\\0&G\end{pmatrix}
\begin{pmatrix}\alpha\\\alpha p\end{pmatrix}.
\]

Here H is the coupled orbital–structure Hessian action, G the nonorthogonal
coordinate metric, and p the microstep. The complete residual remains

\[
r=\alpha(g+Hp-\omega Gp).
\]

We measure r, the predicted gradient, and the true keyframe gradient correction
with their full Euclidean coordinate norms. This is not a norm of their
projection into the AH trial space: components outside that space are included.
The residual zero is unchanged; no inverse metric is needed to calculate it.

For a positive definite G, the norms are equivalent but not numerically equal:

\[
\frac{\|r\|_2}{\sqrt{\lambda_{\max}(G)}}
\leq\|r\|_{G^{-1}}
\leq\frac{\|r\|_2}{\sqrt{\lambda_{\min}(G)}}.
\]

Therefore this is an intentional change of inner progress/acceptance policy,
not an algebraically identical optimization of a physical-norm certificate.
Ill-conditioning may change microstep decisions and the nonlinear trajectory.
The final external energy and orbital-gradient criteria are unchanged and
are evaluated using real accepted-point quantities. No claim of identical
trajectories or universal speedup follows from norm equivalence alone.

The generalized AH equation, root selection (including negative curvature),
HVP, analytic-diagonal preconditioner, orbital chart/retraction, physical step
length, true-energy macro acceptance, and state-average solver are unchanged.
Only the single-state CIAH norm-control path changes. Optional connection
diagnostics and actual preconditioner operations may still use metric inverses;
the change is not a claim that all physical inverses have disappeared.

## Implementation reference

- PySCF [newton_casscf.update_orb_ci](https://pyscf.org/_modules/pyscf/mcscf/newton_casscf.html)
  uses Euclidean norms for predicted and keyframe gradients and their difference.
- PySCF [ciah.davidson_cc](https://pyscf.org/_modules/pyscf/soscf/ciah.html)
  uses the norm of the complete AH residual. It retains the small projected
  overlap matrix; this is not a full-space inverse for a gradient norm.
- Historical commit c2776c4, optimization/neo/solver.cpp, also evaluates the
  complete residual with `.norm()` while retaining the metric in the equation.
  We adopt that separation, not its old FeCl2 endpoint or entire optimizer.

Removed code: deferred physical certification, propagated inverse/defect state,
keyframe triangle-bound logic, and the dense structure-metric audit. No switch
or fallback preserves the old inner inverse path. Unit tests retain generalized
AH dense references and negative-curvature cases; a new nonidentity-metric test
checks a residual outside the trial space and zero physical-norm callbacks,
including after a refreshed gradient.

## Validation

Hanhai22 build/unit job: 1250537. Native parallel array: 1250538, indices
1, 4, 9, 11, 12 (241, FeCl2, QOCGUM, CERRAS, LOFLEA).
32 OpenMP threads, BLAS one thread, exact Libcint integrals, Davidson,
analytic Hessian-diagonal preconditioner, gradient tolerance 1e-3 and energy
tolerance 1e-7 Eh. Frozen binary and source patch accompany the remote build.
No old binary is rerun. Local unrelated RI edits are not transferred.

Compare against saved ciah-bounds-1250459 on Hanhai22 and c2776c4 in
results/history/paired_runs.csv on Hanhai25. The historical FeCl2 endpoint
is incorrect (about 6.825 mEh too high), so its timing is not a correct-energy
speed target. Historical timings are cross-cluster and macrostep accounting
differs; report raw numbers without claiming controlled hardware speedups.

Reproduce the comparison table from saved output (no chemistry calculations):

```sh
python3 results/diagnostics/hanhai22-shared-response-20261012/compare_lean.py \
  ciah-coordinate-1250538 --previous ciah-bounds-1250459 \
  --cases FeCl2 241 QOCGUM CERRAS LOFLEA
```

All five cases converged. The largest endpoint difference from the preceding
correct-energy version is 5.736e-9 Eh. Time is application-reported end-to-end
wall time; RSS is GNU time's process maximum, converted from KiB to MiB.
CPU model is AMD EPYC 9654 for both Hanhai22 series, but node assignments differ.
The input and basis hashes match for every case. These are single measurements,
not a repeated same-node timing study.

| Case | Steps, before → after | HVP, before → after | Time / s, before → after | Peak RSS / MiB, before → after | Final energy / Eh |
|---|---:|---:|---:|---:|---:|
| FeCl2 | 12 → 12 | 232 → 254 | 7.663 → 8.963 | 513.2 → 549.9 | -2181.624464867638 |
| 241 | 3 → 3 | 21 → 21 | 2.301 → 2.272 | 822.0 → 792.4 | -230.720590392829 |
| QOCGUM | 3 → 4 | 37 → 52 | 29.402 → 15.230 | 1969.1 → 1922.5 | -480.094766681858 |
| CERRAS | 5 → 5 | 41 → 51 | 173.956 → 34.043 | 1996.5 → 2026.7 | -397.093435004110 |
| LOFLEA | 4 → 4 | 34 → 34 | 464.043 → 81.751 | 2728.5 → 2696.2 | -399.106646001529 |

Five unit executables passed on Slurm: CIAH, coupled structure, orthogonal
direct-CI, sparse orbital quotient, and Hessian diagonal. Initial build job
1250530 stopped at an obsolete assertion requiring four physical-norm calls;
the assertion now requires zero. Its dependent array 1250532 was cancelled
without running molecular calculations. Job 1250537 passed all tests.

A separate new-binary QOCGUM scalar audit (1250554_15; not a native timing
comparison) records zero `ciah_structure_norm` and `ciah_coupled_norm` events.
The saved old audit recorded 75 structure norm evaluations, 13,787 metric
actions, 13,699 CG iterations, and 18.351 s in those solves. This is eliminated
work, not a faster implementation of the same repeated inverse operation.

The observed improvements are 1.93x for QOCGUM, 5.11x for CERRAS and 5.68x for
LOFLEA. FeCl2 is 17% slower and uses 22 more HVPs. Changing inner norms can
change model-control decisions; this result must not be described as a uniform
performance improvement. No thresholds were tuned in response to these runs.

Historical c2776c4 timings (Hanhai25) are respectively 5.285, 0.916, 6.330,
28.306, and 57.215 s. Thus the present version has not surpassed that baseline
overall. In particular, AH action counts remain above historical values.
The separate FeCl2 historical energy error precludes treating its old runtime
as a numerically equivalent correct-energy benchmark.

Machine-readable values, including historical energy/steps/HVP/RSS, are in
[ciah-coordinate-comparison.csv](ciah-coordinate-comparison.csv). Raw outputs
remain under `ciah-coordinate-1250538/`; all jobs are finished.
