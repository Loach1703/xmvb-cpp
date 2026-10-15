# FeCl2 orbital preconditioner: sampled curvature evidence

## Diagnostic, not a new optimizer

1249034 adds a bounded diagnostic to the retained coherent coarse inverse.
It uses cached orbital directions V and orbital images AV. P is the existing
unshifted analytic-diagonal inverse. No new HVP or integral traversal occurs.
The diagnostic applies P to cached AV and solves a small generalized spectrum.

Define E=V^T A V and T=(AV)^T P (AV). For positive E, solve T u=theta E u.
An exact inverse gives T=E and theta=1. If the full A is positive, these are
Ritz values of A^(1/2) P A^(1/2) in A^(1/2)V. Positivity outside the sampled
space is unproven: this is not a complete-space condition-number certificate.

Sample dimension is at most the existing AH action budget. Shared serial
diagnostic arrays are not replicated per thread. Production scalar tracing
does not run this diagnostic. Operator, candidates and tolerances are unchanged.
Independent tests verify unit spectrum for an exact inverse and correct
reporting of indefinite sampled curvature. Five unit executables passed.
The full-audit run retained 13 macros / 257 HVPs, E=-2181.624464865103 Eh;
diagnostic time is not a native speed benchmark.

| Macro | Sample dimension | theta minimum | theta maximum | Sampled spread |
|---|---:|---:|---:|---:|
| 9 | 21 | 0.4632 | 56.5961 | 122.2 |
| 10 | 22 | 0.5485 | 17.4061 | 31.7 |
| 11 | 22 | 0.4072 | 23.7951 | 58.4 |
| 12 | 25 | 0.3131 | 29.8808 | 95.4 |
| 13 | 28 | 0.1588 | 46.5105 | 292.9 |

The base diagonal does not make all observed curvatures uniform. This is not
evidence of an algebraically incorrect inverse; independent frame tests pass.
Macro 3 has negative sampled curvature and is excluded from the positive
spectral comparison. Negative curvature in the AH operator was not removed.

## Rejected block-coarse experiment

1249037 retained full coupled AH/HVP but calibrated its coarse inverse to
diag(A,C), discarding only Bp/B^Tq in the preconditioner samples. No new HVP,
energy/gradient change, coordinate change, or threshold change was made.
A dense generalized-AH test checked the true coupled operator independently.

| Native run | Macros | HVP | Energy / Eh | GNU wall / s | RSS / MiB |
|---|---:|---:|---:|---:|---:|
| Archived retained 1248947 | 13 | 257 | -2181.624464861511 | 10.44 | 520.6 |
| Block coarse 1249037 | 18 | 389 | -2181.624464822627 | 12.90 | 525.9 |

The ablation is worse. Its change/test were removed locally and from the
execution mirror; no block-coarse switch remains. The generic diagnostic stays
because it supplies evidence without changing production policy.
Orbital-dominated keyframe error does not imply dispensable coupling: orbital
error can arise through structural feedback. Retain the demonstrated coherent
inverse and investigate its complementary-space model and curvature reuse.
Do not simply disable coupling or tune budgets.

Raw logs, native time, hashes, builds and tests are retained in
`fecl2-spectrum-1249034/` and `fecl2-ritz-cost-1249037/`. The archived baseline
was not rerun. Cross-node single-run times are exploratory; increased action
count itself establishes regression.
