# Same-point NEO curvature audit

## Scope and protocol

Hanhai25 compute job `281648` examines NEGPIA at the beginning of point 4,
before the first radius rejection. Accepted energy is -436.257550736293 Eh.
All operator images and finite differences use the same source chart and
the same first core-only candidate direction, normalized in reduced Euclidean
coordinates. Displacements are affine raw-coefficient lifts of that direction.
Sparse gauge canonicalization is disabled for the diagnostic trials, so both
gradient covectors are pulled back through the same immutable source chart.

The diagnostic uses exact Libcint, Davidson and 32 threads. It does not copy
the complete input/integrals: the objective temporarily replaces only orbital
data and restores it. Fixed-structure evaluation uses the accepted structure
coefficients, renormalized in the displaced structure metric. Relaxed evaluation
rediagonalizes the same selected-state problem. Four central-difference steps
check discretization and numerical noise. No full orbital Hessian is assembled.

The dedicated executable exits intentionally with `neo_directional_audit_complete`
after emitting `NEO_COMPONENT_COMPLETE`; application exit 2 is expected and is
not a convergence failure or a converged optimization. The preceding native
baseline in the same job converges. Job `281646` failed to apply the initial
diagnostic patch before any numerical run; `281647` was the first diagnostic,
and `281648` adds the explicit gauge/horizontal checks.

## Decisive measurements

All Rayleigh values below use the same unit reduced-coordinate direction.

| Operator | Directional curvature |
| --- | ---: |
| Core-only | -0.1471097865 |
| Core plus local-active only | -1.4897244781 |
| Complete coupled orbital action | +0.4606028093 |
| Horizontal structure-relaxation correction | -0.0484405931 |
| Complete relaxed HVP | +0.4121622161 |
| Fixed-structure gradient FD, step 1e-4 | +0.4606028153 |
| Relaxed gradient FD, step 1e-4 | +0.4121609730 |

Errors are vector Euclidean errors divided by the larger of one and the
reference vector norm. These are not physical dual-norm residuals.

| FD step | Complete coupled orbital vs fixed gradient FD | Complete relaxed vs relaxed gradient FD | Local-only vs fixed gradient FD |
| --- | ---: | ---: | ---: |
| 1e-3 | 7.05e-7 | 9.74e-6 | 3.43927 |
| 3e-4 | 6.34e-8 | 9.70e-6 | 3.43927 |
| 1e-4 | 8.30e-9 | 1.06e-5 | 3.43927 |
| 3e-5 | 1.54e-8 | 1.46e-5 | 3.43927 |

The fixed-structure error decreases quadratically before numerical noise.
The relaxed check has a roughly 1e-5 floor under production response accuracy;
it does not establish roundoff-level derivative agreement. Independently solving
the horizontal response and adding its adjoint reproduces the relaxed HVP to
1.25e-5 relative error; response relative residual is 5.32e-5. Energy second
differences become noisy at small steps and are not the primary certificate.

## Why the prior local-only interpretation was wrong

`apply_orbital_coupling` enables an additional code path through
`prepare_equal_weight_generalized_eigen_coupling`. This constructs a gauge
coefficient response and selected-space multiplier response, then adds their
orbital adjoint. It enforces the changing structure overlap normalization.
Plain `apply_reduced` with local-active enabled and structure-response disabled
does not take that path.

Denote the core, computational local-active and gauge contributions by
$A_0$, $L$ and $G$, respectively. The correct fixed-structure orbital operator
and relaxed operator are

$$
\bar A=A_0+L+G,\qquad
H_{\mathrm{rel}}=\bar A-B^T C^\dagger B.
$$

Here $B$ and $C$ are the horizontal coupling and structure tangent operators.
The unrestricted computational difference $H_{\mathrm{rel}}-(A_0+L)$
contains both $G$ and the horizontal Schur correction. It must not be called
the horizontal Schur term alone.

On this direction the image norms are: local 5.036, gauge 5.595, horizontal
Schur 0.274, and complete relaxed 1.399. The local-only ablation therefore
omitted a large normalization term. Its poor optimization does not prove that
the physically complete fixed-structure Hessian is a poor approximation.

The production core-only image differs from the relaxed-gradient FD by about
0.992 relative error on this candidate. Solving that approximate subproblem
more accurately cannot remove this model error. The missing terms also change
the directional curvature sign. These are statements about this source point
and direction, not the complete spectrum or every system.

## Next experiment

Use the complete coupled orbital image as an orbital-only approximation,
without solving the horizontal structure-response equation. The existing
`apply_orbital_coupling` API already supplies this image; it also computes
coupling forcing, so this test is not free of overlap/H-S directional work.
It avoids the expensive inverse structure response and keeps the normalization
term. First measure the same six systems before implementing a new general
admission policy or specializing away unused forcing work.

No production numerical code has changed. `neo-component-audit.patch` is a
bounded experiment, and the source mirror and executable were restored.
Raw evidence, input/binary hashes and the exact compiled patch are saved under
`results/studies/locality-20261011/neo-admission-281648/`, with matching files
under `/home/guqqgroup/taoxia/xmvb-runs/native-20261011/` on Hanhai25.
