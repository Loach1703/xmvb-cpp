# FeCl2: residual enrichment without additional operator actions

## Question and fixed controls

The tail's budget exits do not identify an inaccurate small eigensolve.
Native 1249596 has a maximum relative Galerkin residual of 4.03e-12,
whereas the final macro-11 complete residual is 3.65e-4. Its orbital
component is 3.64e-4, and its structure component is 2.91e-5.
The twelve-macro endpoint remains 1.18e-8 Eh above the PySCF coupled
Newton reference. This distinction matters: a projected equation may be
accurate while enrichment of its complement remains slow.

Hanhai22 job 1249619 tests a scalar diagnostic, not a new optimizer.
Coordinates, HVP, preconditioners, keyframe logic, budgets and numerical
tolerances remain unchanged. No new iron convergence task is submitted.
The old native baselines are read, not rerun.

## Geometry of each actually admitted direction

Let r be the pending AH covector residual, omega its Ritz root, G the
physical metric, and q the new direction after existing two-pass Euclidean
orthogonalization. Record

$$
c=\frac{r^Tq}{\|r\|_{G^{-1}}\sqrt{q^TGq}},\qquad
\rho=\frac{q^T(H-\omega G)q}{q^TGq}.
$$

The signed c is the physical angle between the Riesz residual G^{-1}r
and q. It is not the Euclidean angle between two covectors. Its magnitude
is at most one in exact arithmetic; a small value alone does not prove a
poor Newton direction, because a good inverse can emphasize soft modes.
The signed rho supplies that missing directional-curvature context.
Neither number is a whole-space eigenvalue or convergence certificate.

The existing residual norm is cached when r is generated. Gq and Hq are
already evaluated exactly once after admission. The diagnostic therefore
adds scalar dot products only: no HVP, metric action, inverse solve,
retained tall matrix or per-thread carrier. Native scalar diagnostics
report these quantities without changing the direction. Unit counters
check that a residual expansion still makes one HVP/metric evaluation and
no extra dual-norm call.

## PySCF comparison

Reference: PySCF 2.14.0
[davidson_cc](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/soscf/ciah.py)
and
[update_orb_ci](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/mcscf/newton_casscf.py).
The original source was also read from the installed cluster reference.
Both keep admitted directions/images through keyframe RHS changes, expand
with shifted-preconditioned residuals, and permit budget-limited macros.
PySCF's orbital/CI coordinates use a Euclidean norm; our nonorthogonal
coordinates require the complete physical dual norm instead. Replacing
that norm with a diagonal approximation would change the geometry, not
just reduce cost.

In native 1249596, all complete orbital norm calls total 0.2914 s and
structure norm calls 0.1063 s; small spectrum/image construction is
0.0230 s. These measured sections are not the dominant cost of 8.30 s.
The Ritz-only norm timer (0.2484 s) overlaps the complete norm timers and
must not be added again. Expansion HVPs cost 2.2669 s. The audit is intended
to explain their count before proposing any change to the algorithm.

## Status

1249619 completed on anode04; all five unit executables and the original-
guess FeCl2 run exited zero. The run converged in 12 macros/249 HVPs at
-2181.624464872068 Eh, GNU wall 7.40 s, GNU peak RSS 523.69 MiB.
The PySCF coupled-Newton energy difference is +1.1665e-8 Eh. This scalar-
diagnostic run is not a numerical-policy improvement; its timing is on a
different node from 1249596 and must not be advertised as a speedup.

| Macro | Stop | Applied microsteps / budget | AH dimension | Complete residual | Orbital residual | Median absolute physical cosine |
|---|---|---:|---:|---:|---:|---:|
| 9 | micro budget | 21/21 | 23 | 3.20436e-3 | 3.19568e-3 | 0.1295 |
| 10 | micro budget | 23/23 | 25 | 1.00973e-3 | 1.00829e-3 | 0.1296 |
| 11 | micro budget | 24/24 | 26 | 2.86069e-4 | 2.85940e-4 | 0.1296 |
| 12 | micro gradient | 25/28 | 27 | 8.69822e-5 | 8.69255e-5 | 0.1295 |

No tail microstep is clipped. Macro 11 has 99.91% orbital contribution
to squared residual norm. Its two keyframe corrections have captured
squared fractions 3.43% and 3.17%. Admitted directions' final physical
Rayleigh quotients in macros 9--12 are 1.319, 1.494, 2.163 and 1.635.
The consistently small physical alignment supports examining complement
enrichment, but does not by itself certify that these directions are bad
Newton steps or estimate an unrepresented eigenvalue. Do not infer a
whole-space spectral condition number from them.

The complete orbital/structure norm costs in this run are 0.2397/0.0892 s;
small spectrum/image work costs 0.0198 s. The final macro exits at the
unchanged micro-gradient target, not the budget, so "always budget exit"
is not literally true. The tail still consumes many expansions before
reaching that target. No improvement in tail convergence is claimed by
adding a timer or a scalar geometry event.

Build-only 1249623 reran all five unit executables with scalar tracing
enabled, all exit zero. The expansion counter assertion was exercised
with the new geometry events active (16 events in test_ciah). Thus the
no-additional-HVP/metric/inverse-norm assertion is tested in its enabled
mode, not just inferred from the disabled branch.
