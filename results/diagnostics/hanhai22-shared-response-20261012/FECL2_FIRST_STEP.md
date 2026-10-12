# FeCl2 first-step correctness audit (2026-10-12)

Scope: current exact-integral production implementation, FULL OEO, CAS(10,8),
quintet, Davidson, 32 cores on Hanhai22. No archived executable was rerun.
The initial audits change only the diagnostic executable. The subsequent
direction-construction repair enables the existing inner lowest-Ritz-root
search; outer stopping criteria and tolerances remain unchanged.

## Jobs

- 1247319: one NEO accepted step, trace, and sector derivative checks at the
  initial and first accepted point. Exit 2 from control.exe is expected for
  the explicitly bounded one-step run, not a convergence failure.
- 1247328_0..3: higher-accuracy derivative diagnostics on active-containing
  sectors. This is a diagnostic accuracy change, not an optimizer change.
- 1247340: derivative and retraction audit along the observed raw displacement.
  The observed displacement includes finite gauge balancing and is not claimed
  to be the exact original NEO tangent.
- 1247351_0..1: current complete HVP on an independently saved negative direction,
  at the initial and first accepted points. The old direction is data only;
  all actions were evaluated by the newly built diagnostic executable.

## Coordinates and derivatives

All 25 inactive orbitals have 111 tangent columns and 25 gauge columns.
All 8 active orbitals have 103 tangent columns and 33 gauge columns.
Every combined local basis has rank 136. Thus the chart spans the full raw
coefficient space modulo the theoretical FULL/OEO gauge at both audited points;
this is stronger than counting a total of 3599 coordinates.

Coupled orbital/structure derivative checks on inactive-active,
inactive-external and active-external probes have scaled errors of order 1e-9.
These are sampled directional tests, not an exhaustive Hessian validation.
The gauge directional derivative maximum is about 1e-8.

The relaxed HVP discrepancy for active-containing probes at response tolerance
1e-3 disappears when the diagnostic response tolerance is reduced to 1e-8:
the independently assembled coupled-elimination comparison is 1e-11 to 1e-10.
This supports iterative response error rather than an omitted analytic term.
Reported scaled errors use the executable's max(1, reference norm) denominator;
they are not necessarily conventional relative errors for small quantities.

Along the observed first displacement, chart derivative adjoint error is
1.08e-19, retraction tangent relative error is 1.53e-12, complete HVP versus
gradient finite-difference scaled error is 6.17e-9, and gradient versus energy
finite-difference absolute error is 3.13e-8. Newly recorded gradients change
between points; the driver completes the accepted gradient before commit.

## Concrete missing curvature in the first projected model

| Point | Complete-HVP curvature | Gradient-FD curvature | Physical Rayleigh quotient |
|---|---:|---:|---:|
| Initial | -0.00151905668951 | -0.00151905662743 | -0.000741333280858 |
| First accepted | -0.00154823186986 | -0.00154823180402 | -0.000755598709320 |

The known direction remains representable and has a small gradient pairing:
the saved packed direction pairs with raw gradients as -5.25e-7 initially and
+2.31e-7 after the first step. These pairings are not normalized in the current
point's physical metric.

The first accepted projected model instead has a positive minimum eigenvalue
0.0725987988, with 52 orbital directions and one retained structure direction.
It returns TrialReady, not certified KKT convergence: full residual 0.0032654
versus target 0.0010783. The nonlinear keyframe predicts gradient norm 0.0032646
and measures 0.0031282, so the local magnitude guard accepts it. Its energy
reduction is 8.02145e-6 Eh versus prediction 8.00829e-6 Eh.

Therefore a real negative direction exists already at the first point, but
the accepted projected model does not represent it. This is a concrete model
coverage deficiency, not proof of an incorrect complete derivative. KKT
stationarity accuracy alone cannot guarantee selecting the lowest-curvature
branch. Do not fix this by adding endpoint stability convergence conditions.

## Projected-space localization and repair

Job 1247378 reevaluated the retained 52 orbital directions with complete Schur
response, using the same canonically balanced initial point as the NEO trace.
The minimum physical curvature changed from 0.072598798799705 to
0.07220731710594: it remained positive. The known negative direction's metric
coverage was 0.48096068541678. Thus completing structure response alone within
this orbital subspace does not expose the negative branch. An earlier audit
using unbalanced input coefficients (1247373) is invalid and excluded.

The backend previously disabled the solver's lowest-root residual expansion.
Enabling it searches the missing orbital/structure curvature directions during
step construction. This is not an endpoint stability test and adds no outer
convergence condition. Merely adding an independent starting vector (1247383)
did not fix the endpoint; that experimental implementation was removed.

| Calculation | Accepted steps | Final total energy (Eh) | Wall time (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|
| Previous NEO endpoint (historical) | 5 | -2181.617643151932 | — | — |
| Repaired lowest-root search, 1247396 | 17 | -2181.624464422983 | 127.49 | 492.67 |
| PySCF reference (historical) | — | -2181.6244647028957 | — | — |
| Molcas reference (historical) | — | -2181.62446480 | — | — |

The repaired NEO/PySCF difference is 2.7991e-7 Eh, accepted as consistent at the
present convergence settings. Existing response-solver unit tests passed.
The price of lowest-root search is substantial; this establishes correctness
for this case, not a universally efficient search policy.

Job 1247406 tested carrying the lowest-curvature vector instead of the previous
displacement to the next point. It finished in 8.81 s but stopped at
-2181.617645472958 Eh, the undesired high-energy branch. This recycling change
was removed. The retained production change is the validated 1247396 policy.
