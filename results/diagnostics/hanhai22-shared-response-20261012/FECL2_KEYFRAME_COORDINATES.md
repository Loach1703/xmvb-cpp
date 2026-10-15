# FeCl2: coordinates at the actual resumed keyframe direction

## Source comparison

The installed PySCF `newton_casscf.update_orb_ci` commits the incremental
rotation and normalized CI vector using `extract_rotation`, clears `dr`,
and updates the gradient with `g_update(u, ci_kf)`. The Davidson H images
remain attached to the macro-origin model. Its orbital gradient is in the
current rotation frame, not the derivative of one summed additive AO step.

VBSCF currently uses a different, intentional parameterization: orbitals
are an affine additive lift in the immutable origin chart, and raw structure
coefficients are the origin coefficients plus a linear tangent lift.
Normalization is performed inside the energy evaluation. Therefore the
keyframe covector must be pulled back into that same origin chart; simply
substituting a newly constructed chart gradient would be inconsistent.

The orbital projection is the transpose of the same tangent basis used by
`expand_retract_input_tangent`. The structure coordinate lift divides by
sqrt(2*w), while the normalized-coefficient residual pullback multiplies by
sqrt(2*w)/norm. This is the expected derivative of the finite raw-coefficient
normalization, not evidence of an erroneous factor of two.

These are code-level observations, not a proof that the complete numerical
keyframe implementation is correct. PySCF's orthogonal rotation composition
must not be copied blindly into the nonorthogonal VB orbital constraints.

## Numerical discrimination

Job **1248679**, Hanhai22, 32 CPUs, FeCl2 original guess, exact integrals,
analytic-diagonal preconditioner, four macros. Build and `test_ciah` precede
the molecular diagnostic. No optimization threshold or candidate is changed.

Once per macro, take the first ready Davidson microstep immediately after
an accepted keyframe. Normalize it with the existing coupled physical metric
and evaluate the exact fixed-C energy/gradient at the keyframe displacement
plus and minus that direction. Probe spacings 1e-3 and 3e-4 are finite-difference
diagnostics, not optimizer tuning parameters.

`ciah_keyframe_direction_audit` reports:

- keyframe gradient dotted with that actual resumed direction;
- central energy difference and its absolute discrepancy;
- central gradient-difference image in the same origin chart;
- origin Hessian image already cached for that candidate, and their difference;
- orbital and structure contributions to the analytic directional derivative.

The four additional keyframe evaluations make audit timing unsuitable as a
production performance measurement. They do not add an HVP or change the
production microstep. Existing explicit candidate-audit switches gate this
diagnostic; no new environment policy parameter is introduced.

If the energy derivative fails, isolate the incorrect orbital/structure
pullback before fixing it. If it agrees but the local Hessian differs from
the origin image, that difference is finite-model curvature drift, not a
gradient-coordinate bug. It requires a derived model/transport treatment,
not tolerance tuning. No porphine run is submitted before improvement is
demonstrated on this small transition-metal test.

## Completed directional audit

1248679 completed, including a passed `test_ciah`; native molecular exit 2
is the deliberate four-macro cap. The checks are at actual resumed directions,
not another fixed-point energy comparison.

| Macro | Analytic energy derivative | FD error, 1e-3 | FD error, 3e-4 | Local H-image norm | Origin/local H-image difference norm |
|---:|---:|---:|---:|---:|---:|
| 1 | -1.9282522e-5 | 2.63e-9 | 7.10e-9 | 0.11342621 | 0.00281615 |
| 2 | -1.2403900e-4 | 6.12e-10 | 1.00e-8 | 0.01158251 | 0.00957352 |
| 3 | -1.6239364e-4 | 7.00e-9 | 7.25e-9 | 0.01421298 | 0.01093458 |
| 4 | -1.7891766e-3 | 2.98e-9 | 1.25e-8 | 0.07065843 | 0.03171874 |

The sampled scalar derivatives agree at approximately 1e-8 absolute accuracy.
The smaller spacing does not consistently improve the energy derivative,
so this is not an all-directions exact-gradient certification. The gradient
difference images are stable between both spacings. Curvature drift is real
on these directions, but freezing macro-origin curvature is intentional in
PySCF; these observations do not justify replacing the HVP formula.

## Specific incremental-displacement correction

One definite PySCF workflow difference remains: keyframe scheduling uses the
norm of `dr`, which `extract_rotation` commits and then resets to zero. Our
interval used the whole accumulated macro displacement even after refreshing
the gradient. The energy model already owns the post-keyframe incremental
displacement; use that vector for scheduling instead of the macro total.
For the affine origin VB chart this is an exact subtraction/reset operation,
not an adoption of orthogonal orbital rotations or a change of tolerance.

New unit checks assert zero incremental displacement after refresh and that
the subsequent interval contains only the new microstep. Short FeCl2 job
**1248684** completed with passed units and the identical four-step audit
protocol. Both molecular exits are the intended iteration cap, not convergence.

| Implementation | Four-step HVP count | Physical orbital gradient at step 4 | Step-4 energy, Eh |
|---|---:|---:|---:|
| Before interval correction, 1248679 | 80 | approximately 0.0062743 | approximately -2181.6191672 |
| Reset incremental interval, 1248684 | 75 | 0.00825465 | -2181.6191158569 |

The HVP reduction is 6.25%, but the endpoint energy is about 51 microhartree
higher and the gradient is larger. This is not demonstrated convergence
improvement. Retain the mathematically correct interval bookkeeping, but do
not present it as a solution to slow convergence. Iron porphine is not submitted.

## Next identified consistency question

The generalized AH equation uses the complete horizontal structure overlap
metric. The residual and keyframe norms instead use `diagonal_dual_norm`,
which is explicitly a Jacobi approximation, unlike PySCF's Euclidean norm
in its orthonormal CI carrier. Thus the gradient pullback and the norm used
to assess it are different questions. Quantify the effect of this approximation
before changing the model or criteria; do not merely alter the trust constant.
No new production metric inverse or dense structure matrix has been added.
