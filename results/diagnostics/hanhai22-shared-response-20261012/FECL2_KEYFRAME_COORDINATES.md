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

Submission status: pending. No numerical certification or improvement is
claimed before inspecting the actual output.
