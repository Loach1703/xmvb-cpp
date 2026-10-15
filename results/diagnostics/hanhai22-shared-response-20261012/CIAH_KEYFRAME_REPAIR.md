# CIAH keyframe/model correction

## Evidence and scope

Original-guess comparison array 1248419 was cancelled before execution.
The two-step porphine audit 1248292 is reused; no old executable is rerun.
The current source was compared with the complete `update_orb_ci`,
`extract_rotation`, `g_update`, `kernel`, and `davidson_cc` functions in the
installed PySCF reference under `porphine-reference-audit-1247977/python`.

Keeping the macro-origin Hessian images after a keyframe is intentional in
PySCF. A large nonlinear/model discrepancy alone does not establish a wrong
HVP, and retaining those images must not itself be described as a bug.

## Corrections

1. PySCF's converged Davidson yield does not append its old residual. Our
   workspace did. A converged root now retains its existing basis and images;
   a changed gradient is reprojected first and can subsequently generate a
   new residual expansion.
2. The refreshed gradient and the original frozen quadratic energy were
   inconsistent macrostep acceptance models. In particular, the additional
   origin-model prefix restoration does not exist in PySCF's coupled CIAH
   workflow. Remove that restoration and the frozen-model trust-ratio veto.
   Keep the nonlinear relaxed-energy descent check and existing half-step
   backtracking, so this is not an unrestricted acceptance of uphill energies.
3. Refresh the energy anchor together with each accepted keyframe gradient.
   Only displacement since that anchor contributes to its local quadratic
   prediction. Retain the origin H*total-step for separate diagnostics; it
   must not be confused with the refreshed model.

The keyframe energy is for the fixed current structure vector; the final
macro energy is after CI relaxation. A negative fixed-C predicted reduction
does not by itself veto a relaxed-energy descent. Prediction and ratio are
diagnostics, not a new stopping condition. Shortened origin-chart trials use
the origin model only for reporting, not for acceptance.

No gradient/energy tolerance, AH tolerance, step-component cap, keyframe
trust constant, preconditioner, HVP formula, orbital chart, or retraction was
tuned. PySCF's orthogonal orbital rotation is not transplanted into VBSCF.
This is a specific workflow repair, not a claim that every finite-coordinate
transport issue has been eliminated or that convergence is already improved.

## Tests

New `test_ciah` checks cover keyframe energy anchoring, displacement reset,
rejected-step rollback, and converged-root reuse with a nonzero residual
before full-space exhaustion. These are compute-node tests, not local runs.

Hanhai22 job **1248504** builds the new executable and runs unit tests before
concurrently testing FeCl2 and MnF2, four macrosteps each, 32 CPUs per task,
analytic-diagonal preconditioning and exact integrals. Inputs are the same
original guesses as archived `ciah-preconditioners-1247579`; only `itmax`
is changed to four for this bounded diagnostic. Exit 2 means the intentional
cap, not convergence. Pending job 1248492 was replaced before execution.

Both cases enable nonlinear macro auditing. This adds fixed-C gradient work;
their timings are diagnostic timings, not a fair comparison with unaudited
production runs. HVP counts and energy trajectories can be compared with
the archived data. No original-guess long comparison is submitted until
the short checks demonstrate an improvement.

At submission the unit tests and molecular results are **pending**, not
reported as passed or improved.
