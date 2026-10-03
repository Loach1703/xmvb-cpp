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
- [x] Map this independently verified soft orbital mode into the C++ chart
  and compare the full relaxed HVP along it with the C++ energy differences.
  If the curvature is correct, inspect NEO subspace discovery and stopping;
  if it is not, repair the HVP/coupling. Do not replace this diagnosis with
  tighter scalar stopping thresholds or further Davidson parameter tuning.
  Job 268808 (Hanhai25, 32 cores) verifies this direction using dense CI:

  | Contribution / check | Curvature, unit physical orbital step |
  | --- | ---: |
  | Fixed-CI orbital HVP | +0.000003457761 |
  | Structure-response correction | -0.000489596871 |
  | Complete relaxed HVP | -0.000486139111 |
  | C++ retraction energy difference, h=0.01 | -0.000486329554 |

  The complete HVP and energy curvature differ by about 0.039%. Structure
  response relative residual is 3.97e-14. Smaller energy-difference steps
  show cancellation against the ~-2181 hartree total energy; do not use them
  as stronger accuracy certificates. The analytic directional derivative
  is only 4.18e-9 for physical direction norm 1, while physical gradient
  dual norm is 0.000964845: the unstable mode is almost orthogonal to the
  gradient. This validates the sign and response contribution on ONE mode,
  not the entire Hessian. Current NEO initial/final stopping checks only
  first-order stationarity, and its model requests no curvature certificate.
- [ ] Add a separate soft/negative-curvature discovery mechanism near apparent
  convergence, and verify that NEO does not accept this endpoint unchanged.
  Distinguish a discovered negative direction from a global PSD certificate;
  solving the gradient-seeded Newton KKT residual alone cannot establish the
  latter. Use the independently mapped mode above as a strict regression.
  A generic zero-gradient coupled-spectrum probe is now implemented in
  `optimization/neo/curvature.*`. It uses the analytic orbital preconditioner
  only as an SPD spectral metric, does not rotate physical orbitals, and
  does not perform a structure-response solve for every direction. Its
  negative witness must pass a Ritz-residual refinement to the curvature
  scale; finding the first negative projection alone gave poor trial steps.

  **Do not enable the experimental stopping guard by default yet.** The
  backend trial was measured on Hanhai25 job 268813 and then removed from
  the production backend. Its reproducible patch and run artifacts remain
  under `results/diagnostics/pyrrole-molcas-20261002/`:

  | Input / starting point | Original NEO | Experimental guard |
  | --- | --- | --- |
  | FeCl2 dense endpoint | 0 steps, -2181.617644858680 | 20 steps, -2181.623929428968; **not converged**, 114.04 s |
  | F2 OEO original guess | 3 steps, -198.689799668759 | 20 steps, -198.760726276711; **not converged**, 1.05 s |
  | 241 | 6 steps, -230.720590387652, 1.50 s | same steps/energy, 8.62 s |

  The FeCl2 trial escapes the false termination region but leaves physical
  gradient norm 0.02898. It is not a successful final convergence. The F2
  divergence from the old trajectory begins only after the old version
  declares convergence; the new energies are lower, not a forward-energy
  regression. Job 268837 checks the low-energy F2 reference guess: both
  versions reach -198.761111239623 in 5 steps, but the guard adds 45 joint
  actions and increases wall time from 0.20 to 0.53 s. For 241, actions grow
  from 9 to 169 without changing energy. This fails the efficiency target.
  Retain independent curvature discovery as a diagnostic while developing
  reuse/cost control and fixing the slow Newton tail; do not claim a completed
  optimizer fix or silently change all stopping costs.
  Job 268842 validates the retained standalone diagnostic and restored
  production backend. Without using the supplied PySCF direction as its
  spectral seed, the generic probe finds a witness in **7 joint actions**.
  Its orbital component has complete relaxed HVP curvature -0.000648316;
  C++ energy differences at h=0.01 give -0.000648402. The discrepancy is
  about 0.013%. Normal F2 low-guess production is restored to 5 steps,
  -198.761111239623 and 0.20 s. Solver and response-solver unit tests passed.
  The experimental FeCl2 trace shows the next ordinary NEO step consumes
  **486 coupled actions (459 orbital, 27 structure)**. Thus inexpensive
  discovery alone does not fix the subsequent microproblem inefficiency.
  The 241 positive-mode exclusion consumes 80 joint actions. These are
  separate remaining targets; do not blame both costs on initial discovery.

## Required next steps

### 2026-10-03: actual-point gradient / complete HVP verification

- [x] Extend `audit_orbital_gradient` with fixed-chart vector finite
  differences, deterministic random directions, and accepted orbital-table
  overrides. Keep the original raw-coefficient chart at both displacements;
  do not subtract gradients expressed in newly constructed charts.
- [x] On Hanhai25 jobs 269039 and 269054, test FeCl2's dense endpoint and
  the actual first and second accepted points saved by job 268914. Five
  directions (gradient and random) give complete-HVP vector relative errors
  from 9.5e-10 to 6.6e-9 at displacement 3e-5. Gradient energy-difference
  absolute errors reach about 6e-9 to 1.2e-7; cancellation in a roughly
  -2181 Eh energy limits relative accuracy near a stationary point. F2 and
  HAO 241 controls give complete-HVP errors of about 5e-10 and 3e-10.
- [x] Independently assemble the SMALL structure tangent block for diagnostic
  Schur elimination, not as a production algorithm. At all three FeCl2
  points its result agrees with the complete relaxed HVP to 1.6e-15--3e-15;
  structure equation residuals are below 8e-15. This validates the tested
  orbital/structure coupling, not every possible trajectory or direction.
- [x] Test pyrrole 3.5 Angstrom's two-state averaged endpoint in two
  directions. With energy/gradient diagnostic tolerances 1e-12/1e-8,
  shrinking the displacement eventually amplifies Davidson noise. With
  1e-14/1e-10, complete-HVP relative errors decrease from 2.90e-5 at
  displacement 1e-3 to 2.94e-7 at 1e-4 and 2.72e-7 at 3e-5. Energy-gradient
  absolute error reaches 2.96e-8. No full analytic derivative failure has
  been demonstrated at this tested point.
- [x] Direct same-point/same-direction precision ablation, jobs 269063 and
  269068: raw accepted orbital arrays are bitwise identical. FeCl2 default
  Davidson energy differs from dense by only 8.54e-8 Eh, but the raw packed
  gradient L2 error is 5.64e-4 (reference norm about 4.03e-4), and analytic
  HVP vector relative error is 2.89e-4. Diagnostic Davidson tolerances
  1e-12/1e-8 reduce these to 6.88e-8 and 3.98e-9, respectively. These direct
  analytic comparisons are distinct from the much noisier finite differences
  of approximate Davidson gradients. Do not conflate them.
- [ ] Fix the production CI-to-orbital-derivative accuracy contract before
  further optimizer tuning. An energy certificate alone does not certify
  a gradient near small gaps or small orbital gradients. This is evidence
  for selective sensitivity, not a universal RDM formula error. It also does
  NOT explain away the independently established dense-path optimization
  failure or the negative curvature at its old endpoint.
  - [x] Add an independent first-order CI stationarity gate to Davidson,
    measured in Hamiltonian units as residual norm divided by coefficient
    norm. The VBSCF boundary supplies the smaller existing energy/gradient
    budget, instead of inventing a molecule-dependent threshold. Retain the
    relative backward-error and Ritz-energy checks as additional conditions.
    Adding an energy offset cannot weaken this new absolute gate. Expose its
    achieved residual separately from the relative residual in the output.
  - [x] Hanhai25 job 269075: all Davidson unit tests pass, including a
    state with a negligible energy error but a detectable derivative error
    at Hamiltonian offsets 0 and -2000. At the exact same FeCl2 raw orbitals
    and packed direction as the dense reference, the default gradient L2
    error drops from 5.6436e-4 to 2.2141e-7; full analytic HVP relative error
    drops from 2.8913e-4 to 2.8345e-6. The energy is
    -2181.6176448586803 Eh. No orbital derivative formulas were changed.
  - [x] Same-node, same-input, same-backend comparison against HEAD 94d343d,
    Hanhai25 job 269080 (32 threads, default Davidson, exact integrals):

    | Input | Optimizer | Steps before / after | Wall seconds before / after |
    | --- | --- | ---: | ---: |
    | F2 low guess | block-LBFGS | 15 / 15 | 0.108 / 0.113 |
    | F2 low guess | TNHVP | 8 / 8 | 0.250 / 0.246 |
    | F2 low guess | NEO | 5 / 5 | 0.198 / 0.212 |
    | HAO 241 | block-LBFGS | 9 / 9 | 0.991 / 1.188 |
    | HAO 241 | TNHVP | 6 / 6 | 1.612 / 1.783 |
    | HAO 241 | NEO | 6 / 6 | 1.497 / 1.706 |

    All six pairs converge with unchanged outer step counts. Final energies
    change by at most 2.2e-9 Eh. This is an accuracy repair, NOT a speedup:
    241's additional CI work costs roughly 11--20% in these single runs.
    The initial job 269075 used input-default TNHVP; do not compare its
    8-step F2 result with a prior 5-step NEO run as an algorithm regression.
    Final-source validation job 269085 passes Davidson, NEO, response-NEO,
    and eigen-response tests and reproduces the same FeCl2 derivative errors.
  - [x] Test final optimization energy, not just same-point derivatives:
    Hanhai25 job 269147_1 uses revision 35d0871, the original FeCl2 guess,
    exact integrals, default Davidson and block-LBFGS. It reports convergence
    after 256 steps, 67.12 s, at -2181.6176448746 Eh. Physical gradient dual
    norm is 9.88834e-4. The converged same-guess PySCF Newton reference is
    -2181.624464883727796 Eh (job 249623, 24 macro iterations), leaving
    a gap of **0.0068200091 Eh**. Thus the CI stationarity repair does NOT
    resolve the main energy discrepancy. Do not label this energy validated.
  - [x] Companion NEO job 269147_0 was stopped after 1312 steps at
    -2181.6176414357 Eh. This is an intermediate value, not convergence.
    Both jobs retain the original input's 2000-step budget, not a 30-step
    diagnostic cutoff. Logs are under
    `results/diagnostics/pyrrole-molcas-20261002/ci-convergence/`.
  - [ ] A CI stationarity residual is NOT a general bound on the orbital
    derivative error: inverse structure gaps and orbital/structure coupling
    still matter. The new gate repairs the demonstrated premature acceptance,
    not every possible near-degenerate case. A complete sensitivity-aware
    derivative certificate remains separate work; do not mark it complete
    on the strength of this single-point regression.

Logs and the complete four-displacement error grid are under
`results/diagnostics/pyrrole-molcas-20261002/vector-audit/`;
`derivative_errors.csv` is generated by `summarize_vector_audit.py`.
The initial diagnostic's `vector_fixed_ci_relative_error` and
`vector_response_relative_error` fields must NOT be used as component-error
certificates: holding normalized structure coefficients fixed in a Rayleigh
quotient is not the same split as the internal fixed-CI/response Hessian.
That incorrect diagnostic comparison was removed and replaced by the explicit
small-block Schur reference. No production derivative formulas were changed.

### 2026-10-03: coupled-response certificate audit

- [x] Fix false global-curvature certification: a complete orbital basis and
  a converged lowest Ritz residual are insufficient when the structure basis
  is incomplete. The regression has a retained curvature of -1 and an omitted
  Schur curvature of -8. Both bases must be complete before reporting the
  stronger global certificate.
- [x] Add candidate-local structure-residual correction before accepting a
  coupled KKT certificate. Solve `C dq = -r_s`, evaluate `B^T dq` and `C dq`,
  and check both corrected residuals. A missing response direction is added
  to the retained structure basis only when needed. Do not solve response for
  every retained orbital direction. The regression uses `C = 1e-6`: an
  apparently acceptable `|r_s| = 5e-5` propagates into an orbital error of
  0.05; the corrected solution agrees with the explicit Schur reference.
- [x] Validate the production correction on Hanhai25 job 268983 (32 threads,
  one-electron preconditioner, exact integrals). NEO and response-NEO unit
  tests pass. F2 low-energy reference remains 5 steps at
  -198.761111239623; 241 remains 6 steps at -230.720590387652 (dense) and
  -230.720590385933 (Davidson). Dense wall time changes from 0.221 to
  0.251 s for F2 and 1.603 to 1.678 s for 241. Davidson wall time changes
  from 0.216 to 0.234 s and 1.278 to 1.552 s, respectively. These are
  single-run checks, not a speedup claim. Additional response work is counted.
- [x] Reject the simultaneous orbital/structure expansion trial: on the
  identical FeCl2 second-step fixture, actions change from 459 orbital + 27
  structure to 413 + 84. It does not improve total action count. Canonical
  coordinate completion is not responsible for the original 486 actions.
- [ ] Resolve the remaining FeCl2 optimizer failure. Analytic diagonal alone
  reduces the diagnostic second-step action count from 486 to 104, but its
  full trajectory still times out after 29 accepted steps in 300 s. With
  the new candidate-local response certificate, the experimental negative-
  curvature guard plus analytic diagonal still times out after 6 steps in
  180 s (last printed energy -2181.6240866197). Neither is a successful
  convergence. The guard/diagonal default change is NOT retained in the
  production backend; its patch and logs remain under results/diagnostics.
- [ ] Derive metric-consistent residual/error control and investigate the
  accepted-point coordinate conditioning. Do not blindly remove the absolute
  floor or repeat system-specific tuning. Small coordinate-space structure
  residuals do not bound the eliminated orbital error when C is ill-conditioned.
  The production correction addresses that certificate defect, not the entire
  nonlinear convergence problem.
  **Confirmed long-tail mechanism, not yet an optimizer repair:** baseline
  revision 94d343d also stagnates at step 642, -2181.6176426469 Eh
  (job 269150, diagnostic wall-time limit). Do not blame the entire failure
  on the CI stationarity gate; retain that validated accuracy change.
  Job 269164 streams 689 accepted-step records with no change to the NEO
  algorithm. At step 280:

  | Quantity | Value |
  | --- | ---: |
  | Source coordinate gradient L2 | 2.3690e-4 |
  | Source physical gradient dual norm | 1.5495e-3 |
  | Accepted coordinate KKT residual | 2.6822e-4 |
  | Coordinate KKT target | 1.000056e-3 |
  | Physical orbital KKT residual | 1.5542e-3 |
  | Accepted physical gradient dual norm | 1.5556e-3 |
  | Orbital HVP / structure actions | 1 / 42 |
  | Trust radius / shift / boundary | 1 / 0 / false |
  | Actual decrease / trust ratio | 4.0973e-10 / 0.9941 |

  The coordinate KKT gate certifies a step while its physical residual is
  essentially the entire physical gradient. The outer physical-gradient
  gate correctly refuses convergence, and globalization accepts tiny
  positive decreases. The next keyframe repeats one orbital direction and
  dozens of structure actions. At step 689 the physical residual remains
  1.3475e-3 and the decrease is only 1.69e-10 Eh. This is not trust-radius
  exhaustion or evidence of inaccurate CI energies. Structure-response
  correction alone does not repair the missing orbital accuracy.
  The new observer flushes `--neo-trace` per accepted step; physical norm
  diagnostics run only when that observer is enabled. NEO/response-NEO tests
  pass; F2 and 241 remain 5 and 6 steps at their previous energies.
  **Next:** use one physical residual/forcing contract for inner and outer
  decisions and enrich the orbital space when that certificate fails.
  Do not simply remove the scalar floor, tune a smaller fixed threshold, or
  claim this also resolves the independent negative-curvature escape problem.
- [ ] Do not repeat the already established fixed-orbital energy comparison.
  Forward agreement is established; it does not separately prove every RDM
  component or every derivative. A real negative-curvature HVP direction also
  agrees with energy differences. Further tests should target the actual
  failing derivative directions and optimizer certificates, not assume RDM
  failure solely because the system is ill-conditioned.

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
