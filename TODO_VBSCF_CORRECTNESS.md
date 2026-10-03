# VBSCF correctness: FeCl2 and pyrrole

Status: **neither energy discrepancy is resolved**. Do not label these cases
correctly converged from the printed stopping reason alone. Do not change the
wavefunction constraints merely to reproduce a reference energy.

## Verified on Hanhai25, 2026-10-03

- [x] Certify the actual generated FULL spaces, not just their counts.
  `test_full_space` checks spatial-occupation sectors, highest-weight spin
  purity, and exact modular incidence rank. FeCl2 CAS(10,8), quintet:
  308 sectors and 420 independent structures. Pyrrole CAS(11,10), doublet:
  8,350 sectors and 27,720 independent structures. Job 268767 passed;
  the CMake target also passed in job 268782.
- [x] Compare physical gradients at both exported endpoints using independent
  PySCF gradients and central energy differences, without optimizing again.
  Job 268768: FeCl2 physical flag dual norm 0.0017390076; pyrrole warm-restart
  endpoint 0.0032992474. Reported production values were respectively
  0.0009312280 and 0.0009828765. The import/export and CI solve paths differ;
  this discrepancy alone does not establish which layer caused it.
- [x] Compare C++ derivatives with its own energy differences along the
  independently constructed physical direction. Job 268782:

  | FeCl2 structure solve | Analytic directional derivative | Central difference, physical step 0.0003 |
  | --- | ---: | ---: |
  | dense | -0.0017388357 | -0.0017402832 |
  | default Davidson | -0.0028275337 | -0.0017457122 |
  | reference-accuracy Davidson | -0.0017387477 | -0.0017402832 |

  Dense and reference-accuracy Davidson agree; default Davidson is materially
  inaccurate for this derivative although its energy differs from dense by
  only about 0.000000108 hartree. Finite steps retain truncation error; do not
  interpret every small derivative difference as a code defect.
  Pyrrole reference-accuracy Davidson gives directional derivative
  -0.0033005291 and physical norm 0.0033032490, close to independent PySCF.
  These are diagnostic accuracy settings, not a proposed production cutoff.
- [x] Repair the state-average audit tool: its history replay, trial checks,
  dense comparisons, and displaced gradients must use the same selected
  states and weights as the accepted point, not hard-coded root zero.
  Job 268784 passed a dense F2 two-state audit: directional HVP finite-
  difference relative error 3.95e-9.
- [x] Isolate diagonalization from orbital optimization using the original
  FeCl2 guess and block-LBFGS, changing only the eigensolver. Job 268789:
  dense takes 258 accepted steps / 287.83 s to -2181.617644858690;
  Davidson takes 117 / 31.23 s to -2181.617641690341. The dense result remains
  0.006820025 hartree above the same-guess PySCF Newton reference. Thus
  Davidson derivative errors are real but do not explain the principal gap.
  Both runs use the one-electron block preconditioner.
- [ ] Audit accepted-step gauge changes against block-LBFGS secant history.
  The Armijo routine canonicalizes accepted orbitals and exposes
  `accepted_chart_changed`, but the block backend does not request that flag
  and clears history only when the chart rank changes. FULL OEO disables
  active balancing, not inactive balancing. This is a code-level risk; its
  occurrence and impact in the FeCl2 trajectory have not been measured.
- [x] Probe the actual dense endpoint with an independent matrix-free coupled
  Hessian and relaxed-energy directional differences. Keep the action budget
  bounded; absence of a negative Ritz value is not a positivity certificate,
  and negative curvature at a nonstationary point does not establish a saddle.
  Job 268798: independent orbital gradient L2 is 0.000970395, consistent
  with the reported physical norm 0.000964406. Unscaled 40-action probes
  found positive Ritz values but large residuals; they do not certify
  positive semidefiniteness. Job 268799: PySCF joint Newton restarted from
  exactly this occupied-space endpoint reaches -2181.622437040922 after
  8 macro iterations (17 keyframes, 88 HVPs), decreasing energy by
  0.004792182 hartree while retaining S^2=6. This bounded restart is NOT
  converged. It proves substantial further descent is possible and rules
  out labeling the dense endpoint a validated local minimum.
  Job 268803 adds a diagonal congruence to expose the soft spectrum:
  a 60-action projected joint probe finds a negative Rayleigh direction.
  The scaled Ritz value is -8.451875 (residual 8.4e-9); the unscaled
  normalized joint Rayleigh value is -0.000174981. These are different
  coordinate metrics, not interchangeable Hessian eigenvalues. Independent
  CASCI-relaxed energy differences along its normalized orbital component
  give curvature -0.000972168 at h=0.01 and -0.000972957 at h=0.003.
  Smaller displacements show cancellation, so use the two larger probes
  as the sign check. This establishes real negative local curvature at the
  dense endpoint, not a formally stationary saddle or a global-minimum claim.
- [ ] Map this independently verified soft orbital mode into the C++ chart
  and compare the full relaxed HVP along it with the C++ energy differences.
  If the curvature is correct, inspect NEO subspace discovery and stopping;
  if it is not, repair the HVP/coupling. Do not replace this diagnosis with
  tighter scalar stopping thresholds or further Davidson parameter tuning.

## Required next steps

1. [ ] Certify CI error propagation into the physical orbital gradient at the
   in-memory accepted point, including recycled eigenvectors and gauge changes.
   An energy-based CI certificate is not a first-order derivative certificate.
2. [ ] Implement an error-controlled CI/gradient accuracy contract. Derive the
   error estimator before choosing the policy; do not substitute a universally
   tiny Davidson tolerance for that contract.
3. [ ] Verify full relaxed HVP and orbital/structure coupling on the real
   single-state and state-average inputs using identical states at every probe.
4. [ ] Integrate the derivative certificate into convergence decisions and
   Newton/NEO model construction for all relevant optimizers.
5. [ ] Re-run both systems from their original guesses. Report state energies,
   physical stationarity, CI residuals, iterations, actions, time, and RSS.
   If stationary solutions still differ, resolve state/branch selection with
   cross-restarts and curvature evidence rather than asserting equality.

Implementation references: `tests/input/test_full_space.cpp`,
`src/tools/audit_orbital_gradient.cpp`, and `src/tools/audit_newton_step.cpp`.
Detailed input, direction, endpoint, and finite-difference artifacts are under
`results/diagnostics/pyrrole-molcas-20261002/gradient-comparison/`.
