# FeCl2 tail: projected solve versus complementary residual

## What the archived trace establishes

Run 1248805 uses the original guess, exact integrals, 32 bound Hanhai22
cores, the full-frame analytic diagonal, and unchanged CIAH policies.
It converged in 16 macros / 334 orbital HVPs, 11.87 s, 516.8 MiB peak RSS,
at -2181.624464776006 Eh. The archived PySCF reference is
-2181.6244647028957 Eh; their difference is -7.3110e-8 Eh.

Macros 9--16 all stopped at `micro_budget`, not component clipping. Each
had two nonlinear keyframes. Their final predicted gradient norms were
0.003682, 0.001898, 0.002273, 0.001897, 0.001019, 0.000792, 0.000850 and
0.000437. None reached the existing inner gradient stop (0.00009).
The final AH eigenvalue at macro 16 was -7.2475e-9; its AH residual norm
was 0.0004368200504, almost the predicted gradient 0.0004368200655.
The AH shift was therefore not removing that observed residual.

For these tail queries, the Euclidean retained-space residual divided by
the full Euclidean residual was below 2.222e-13. This certifies the
Galerkin solve **inside the retained space**, not convergence of the full
AH problem and not a physical-norm decomposition. More accurate dense
projected eigensolves cannot remove a complementary-space residual.

Native scoped timing counters recorded:

| Operation | Accumulated seconds |
|---|---:|
| Projected H/metric formation and checking | 0.1867 |
| Small spectrum solve and step/image combination | 0.0226 |
| AH residual dual norms | 0.7182 |
| All traced structure dual norms | 1.2021 |

The last two rows **overlap**; do not add them. All structure norms made
36,780 overlap actions in 808 solves, with certified relative residuals
at most 1.4893e-8 (the unchanged metric solver accuracy). These are explicit
wall-clock scope measurements, not CPU-profiler attribution percentages.

## PySCF comparison

Inspected installed PySCF sources on Hanhai22:
`pyscf/soscf/ciah.py::davidson_cc` and
`pyscf/mcscf/newton_casscf.py::update_orb_ci`.
Both retain basis vectors and Hessian images while refreshing the gradient
RHS. After an unconverged yield, PySCF expands the previously yielded
residual; it does not recompute a new residual before each resumed expansion.
Our generalized version uses the physical metric and its dual residual norm,
rather than PySCF's Euclidean norm for orthonormal coordinates.
The gradient-dependent microstep budget is also present in PySCF's coupled
driver. A budget exit is a macro refresh policy, not a certificate that the
Newton equation has been solved. Raising that budget is not the investigated
fix. The previous immediate-keyframe RHS-refresh experiment was already
rejected (21 macros / 390 HVPs versus 18 / 353).

## Rejected reference-frame experiment

1248875 rotated only the diagonal preconditioner's virtual reference frame
to diagonalize the accepted inactive Fock operator. Real nonorthogonal
orbitals, the chart, HVP, retraction, step cap and tolerances were unchanged.
Five unit tests passed, including exact/RI diagonal contractions and the
independent frame inverse. It converged, but did not improve efficiency:

| Run | Macros | HVP | Energy / Eh | GNU time / s | RSS / MiB |
|---|---:|---:|---:|---:|---:|
| 1248805, retained | 16 | 334 | -2181.624464776006 | 11.87 | 516.8 |
| 1248875, reference-frame experiment | 16 | 340 | -2181.624464793568 | 12.05 | 484.1 |

The experiment's production and test changes were removed. Lower RSS in
one native run is not sufficient reason to retain a failed convergence
optimization. No alternate switch or unused canonicalization path remains.

## Next bounded diagnostic

1248881 restores the retained implementation and adds opt-in scalar timings
for residual preconditioning, two-pass orthogonalization, metric images and
HVP images. It also records the fraction of each preconditioned direction
remaining outside the retained Euclidean basis. No additional operator
action, dense Hessian or large thread-local array is introduced.
This separates expensive images from redundant expansion directions before
designing the next preconditioner/subspace change. No new porphine long run
was submitted during this investigation.

The diagnostic finished on anode39 with all five unit tests passing. It
reproduced 16 macros / 334 HVPs and -2181.624464776000 Eh, GNU wall 11.82 s,
RSS 507.3 MiB. No efficiency gain is claimed for instrumentation.

| New scoped measurement | Accumulated seconds |
|---|---:|
| Residual preconditioner | 0.1206 |
| Two-pass Euclidean orthogonalization | 0.0162 |
| Admitted-direction metric action | 0.1803 |
| Admitted-direction Hessian action | 2.6149 |
| Projected formation/checking | 0.1856 |
| Spectrum solve and step/image combination | 0.0221 |
| All structure dual norms | 1.1880 |

These scopes do not account for the entire 11.82 s; accepted-point builds,
nonlinear keyframe evaluation, CI relaxation and other work remain outside
them. HVP scope includes the coupled action, not just orbital arithmetic.
Do not interpret any row as a hardware-profiler CPU percentage.

In macros 9--16, the smallest fraction of a preconditioned direction
remaining after orthogonalization ranged from 0.5231 to 0.7263. Expansion
is not stalled by nearly dependent directions. All eight macros still
ended at the unchanged micro budget, with no component clipping. The
failure is inefficient **full residual contraction despite independent
directions**, rather than an inaccurate projected eigensolve, direction
collapse, or virtual-frame orientation alone. This narrows the next
investigation to the spectral effectiveness of the preconditioned full
operator and the gradient changes introduced at nonlinear keyframes.
It does not yet prove which of those mechanisms is dominant.
