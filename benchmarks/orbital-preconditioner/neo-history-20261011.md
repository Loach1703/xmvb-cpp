# Historical NEO candidate strategy

## What changed since the measured good version

The saved `c2776c4` executable reproduces MnF2's 199 orbital HVPs on
Hanhai25 (job 281440). Its 11 reported macros contain **14 accepted exact-gradient
points**. The current baseline has 15 accepted points and 522 HVPs.
The historical and current trust-radius update formulas are identical.

Three differences must not be confused:

1. Historical coupled solves stopped on coordinate residuals and used the
   coordinate gradient for the nonlinear forcing term. Current certified
   solves use the physical dual norm, including the candidate-local structure
   response check. Removing these checks is not a correctness-preserving rollback.
2. Historical shifted preconditioning solved the local positive model to
   roundoff-level accuracy. Restoring that policy alone on current operators
   gives 624 HVPs, 14 accepted points and 14.817 s in job 281441. It was removed.
   This experiment does not identify shifted preconditioner accuracy as the
   cause of the favorable old trajectory.
3. The current early candidate policy was boundary-only. An interior candidate
   could not be tried until physical KKT certification completed, even if
   the historical coordinate criterion would already have returned it.

Job 281440 started while the precise-preconditioner executable was being
rebuilt. Its historical rows are valid historical reproductions; its updated
rows are **not a uniform baseline panel**. MnF2's updated hash is the precise
preconditioner executable. The independent matched panel below uses immutable
baseline `72ab2e5f5ada70f2af9964c6e86a32b66fbaf83e1cd558b7180fddcaf14117d0`.
The batch script now snapshots the updated executable before running either
variant, so later builds cannot change the executable selected within a pair.

## Recovery scope

Reuse the historical coordinate-residual condition as **candidate readiness**,
not as convergence. It uses the original coordinate-gradient forcing rule
and existing absolute tolerance, introducing no molecule-specific threshold.
It can yield interior as well as boundary candidates. The solver returns
`TrialReady`, never `Converged`, for this condition.

The actual trial energy must pass the existing trust-region acceptance rule;
the completed relaxed gradient must pass the existing physical-norm trust
guard. The accepted point gets a new exact gradient and local model. A failed
gradient guard disables early trials and resumes the retained space at the
same point/radius. Repeated yields require a new independent basis direction.

Retain the existing physical-gradient proximity gate. Additionally, disable
early candidates once the coordinate gradient infinity norm passes the
requested outer threshold. At that point a remaining physical gradient must
be resolved by a certified coupled solve, not by the historical coordinate
criterion. Without this second gate, FeCl2 stagnates near physical gradient
0.0034 with tiny energy decreases (job 281447_2, stopped after diagnosis).

The physical KKT certificate, structure-response certificate, outer energy
threshold, OEO physical-stationarity gate, exact/RI operators, tile storage,
and OpenMP data ownership are unchanged. No new preconditioner, large array,
HVP replay, compatibility backend or user parameter was added.

This is a coordinate-dependent historical **trial heuristic**, not a new
coordinate-invariant inexact-Newton convergence theorem. It does not promise
quadratic convergence of its uncertified intermediate steps.

## Matched native results

Job 281453 built the change and passed `test_neo`,
`test_neo_response_solver` and `test_shifted_metric_preconditioner`.
Job 281454 ran five pairs independently in parallel, each on the same node,
with identical inputs, 32 compact-affinity OpenMP cores, one BLAS thread,
exact Libcint integrals and Davidson. Thresholds remain gradient 1e-3 and
adjacent energy change 1e-7 Eh. These are application end-to-end wall times
and GNU-time process peak RSS, not profiled times.

| System | Accepted points old/new | Orbital HVP old/new | Structure actions old/new | Wall s old/new | RSS MiB old/new |
| --- | ---: | ---: | ---: | ---: | ---: |
| MnF2 | 15 / 14 | 522 / 180 | 1748 / 486 | 9.493 / 4.385 | 505.1 / 493.8 |
| FeCl2 | 3 / 3 | 302 / 287 | 1145 / 1008 | 6.244 / 5.829 | 415.6 / 401.5 |
| 241 | 6 / 6 | 10 / 10 | 119 / 119 | 2.143 / 2.318 | 938.5 / 944.9 |
| 240 | 8 / 8 | 13 / 13 | 0 / 0 | 3.520 / 3.655 | 1743.6 / 1688.6 |
| F2 | 5 / 5 | 11 / 11 | 25 / 20 | 0.143 / 0.136 | 57.9 / 62.0 |

MnF2 uses 65.5% fewer HVPs and 72.2% fewer structure actions; its measured
wall time falls 53.8%. The other small timing/RSS changes are single-run
observations, not significant universal speedups. 240 is orbital-only and
does not exercise the new coupled candidate criterion.

| System | Baseline energy Eh | Recovered-policy energy Eh | Difference Eh | Final physical gradient old/new |
| --- | ---: | ---: | ---: | ---: |
| MnF2 | -1348.893353132275 | -1348.893351870672 | +1.261603e-6 | 8.4829e-4 / 6.1188e-4 |
| FeCl2 | -2181.617644857094 | -2181.617645215967 | -3.58873e-7 | 9.1958e-4 / 9.8856e-4 |
| 241 | -230.720590376551 | -230.720590376550 | +1e-12 | 4.5057e-4 / 4.5057e-4 |
| 240 | -343.446302348804 | -343.446302348805 | -1e-12 | 1.4854e-4 / 1.4854e-4 |
| F2 | -198.751155829774 | -198.751155829492 | +2.82e-10 | 1.1258e-4 / 1.2817e-4 |

**Do not claim equal-energy acceleration from this panel.** Both variants
satisfy the requested stopping tests, but MnF2 and FeCl2 endpoints differ by
more than 1e-7 Eh. The 1e-7 adjacent-step test is not an absolute energy-error
certificate. The difference is recorded rather than dismissed as harmless
or hidden by the speedup.

An additional, explicitly stricter diagnostic (job 281456, gradient 1e-4,
energy change 1e-8) does not yet establish equal endpoints. Its MnF2 baseline
fails after 16 accepted points in the projected generalized-eigen response
MINRES certificate; the batch consequently does not run the updated variant.
This is a failure of the unchanged baseline at the stricter setting, not
evidence for or against the new candidate policy. The FeCl2 baseline also
fails in the response MINRES certificate, after reaching -2181.623943989361 Eh;
it is not a converged reference energy. Neither updated variant runs after
the baseline failure. These diagnostics are not included in the converged
native comparison above and do not certify identical endpoint precision.

The regression fixtures check boundary-trial continuation with cached HVPs
and interior historical trial readiness on an anisotropic metric. An interior
trial fails the physical certificate, while strict continuation recovers the
explicit Newton solution without repeating retained HVPs. Job 281462 passed
all three tests after rebuilding the final fixture separating CIAH's
three-direction test from the historical interior test. Its native MnF2
run reproduces 14 accepted points and 180 HVPs.

## Evidence

Local raw data: `results/studies/locality-20261011/neo-early-281440/`,
`neo-expansion-281441/`, `neo-expansion-281448/`, `neo-early-281454/` and
`neo-early-281456/`. The parser `neo_history_results.py` extracts energies,
accepted-point counts, actions, times, RSS and input/binary hashes directly
from saved output and TSV files. `history-matched.csv` contains job 281454's
ten converged runs. No workload was executed locally or on a login node.
