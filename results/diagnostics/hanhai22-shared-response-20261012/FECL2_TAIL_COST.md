# FeCl2 CIAH tail: subspace convergence versus microstep exhaustion

## Scope and provenance

Hanhai22, Slurm 1249441, anode31, 32 bound cores; original FeCl2 input,
exact integrals, Davidson, analytic Hessian diagonal. Numerical controls and
operators are unchanged from retained commit ee8d2dd5. This revision adds
scalar timing of existing forward phases and diagonal construction only.
Archived binaries were not rerun. Five unit executables passed.

Native output reports convergence in 13 macros / 254 orbital HVPs,
energy -2181.624464847672 Eh, GNU wall 9.02 s, peak process RSS 523.4 MiB.
The endpoint is 3.6061e-8 Eh above the PySCF coupled-Newton reference
-2181.624464883733 Eh. The different action count from the archived 14/280
run is not attributed to timing instrumentation or claimed as an algorithmic
improvement; no new numerical policy was implemented here.

## Measured costs

These are instrumented native wall scopes, not sampling percentages. They
are not an exhaustive partition of total wall time. Nested scopes must not
be added indiscriminately: structure norm time is included in Ritz norm time,
and HVP expansion includes its image-observation/coarse-space ingestion.

| Scope | Accumulated seconds |
|---|---:|
| Coupled HVP expansion, including image observation | 2.12784 |
| Analytic orbital diagonal construction, 14 builds | 1.79730 |
| AH residual preconditioning | 0.17191 |
| AH projection and metric-spectrum checks | 0.12687 |
| AH augmented spectrum and image reconstruction | 0.01443 |
| AH residual dual norm | 0.22477 |
| All complete structure dual norms (nested subset) | 0.10919 |
| Orbital preparation, 36 forward/keyframe evaluations | 0.14635 |
| AO effective one-electron forward | 0.09966 |
| Active one-electron forward | 0.00036 |
| Active two-electron forward | 0.20323 |
| Structure forward setup | 0.07186 |
| Forward CI eigensolver | 0.11027 |

Structure Riesz solves use at most one CG iteration, 856 overlap actions,
and maximum certified relative residual 2.0083e-15. The complete norm remains
certified against the true metric; no precision was relaxed.

VTune software sampling was attempted in a separate bounded eight-macro
capture, but collection was denied by ptrace_scope. Its application did not
run, so there are no usable hotspot percentages. The batch completed but
the profiling step failed; retain profile/run.err rather than treating it as
a successful capture.

## Why the tail exits at the budget

The current PySCF coupled-Newton implementation has the same ordering:
an AH candidate is accepted as a usable approximate microstep, the predicted
gradient is updated, then the gradient-dependent microstep budget is checked
before deciding whether another keyframe is needed. Its budget is a model
refresh policy, not a certificate of exact Newton-equation convergence.
It also retains prior Hessian images and uses the yielded residual for the
next expansion. Appending that residual after a keyframe is not, by itself,
a discrepancy from PySCF.

Our projection is solved accurately *within its current basis*: the largest
Galerkin residual divided by the full residual is 2.47e-12. This does not
mean the full AH residual is small. The latter is mostly outside the trial
subspace and still requires effective new orbital directions.

| Macro | Exit | Full last AH residual | Orbital share of residual squared | Clipped microsteps |
|---|---|---:|---:|---:|
| 11 | micro_budget | 1.93192e-3 | 98.70% | 0 |
| 12 | micro_budget | 2.54769e-3 | 99.47% | 0 |
| 13 | micro_budget | 4.20293e-4 | 99.34% | 0 |

The final budget is 25 microsteps, with a 27-vector AH space. The last full
model gradient, 4.20293e-4, remains above the unchanged inner target 9e-5.
The projected eigensolver is therefore not falsely declaring a full-space
Newton solution. These tail steps are not limited by the component clip.

Keyframe mismatch is also material: maximum mismatch/predicted norm is
0.88, 2.92 and 0.30 in macros 11, 12 and 13. Approximately 99% of the
mismatch squared in these keyframes is orbital. A tiny AH eigensolve or an
already-fast metric inverse cannot fix that model discrepancy.

Supported next target: orbital expansion/preconditioning and the frozen
orbital model across true-gradient keyframes. Diagnose these separately;
do not increase the budget, loosen residuals, or interpret Galerkin accuracy
as a full-space residual certificate. No new iron-porphine long run was
submitted on the basis of this instrumentation-only result.

## Rejected trial

1249431 replaced the raw structure AH diagonal by a mean normalized
Rayleigh gap times the exact horizontal structure metric inverse. It
converged in 15 macros / 271 HVPs, -2181.624464860973 Eh, 9.62 s,
531.3 MiB. The retained archived comparison was 14/280, 9.81 s on a
different node. There is no demonstrated overall improvement; all trial
production and unit-test code was deleted, with no option or fallback.

## Sources and reproduction

- [PySCF coupled-Newton update_orb_ci](https://github.com/pyscf/pyscf/blob/master/pyscf/mcscf/newton_casscf.py).
- [PySCF davidson_cc](https://github.com/pyscf/pyscf/blob/master/pyscf/soscf/ciah.py).
- Native evidence: fecl2-ritz-cost-1249441/native/{output.out,run.err,neo.tsv,resource.time,hashes.txt}.
- Programmatic extraction: summarize_ritz_cost.py RUN_DIRECTORY.

The source links refer to the implementation inspected on 2026-10-13;
these moving master links are not a pinned PySCF benchmark version.
