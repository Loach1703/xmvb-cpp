# Porphine CIAH correctness audit — 2026-10-12

Scope: spin space, derivatives, and macro/keyframe consistency. No changes to
optimizer thresholds, no performance optimization, and no repeated comparison
of fixed-point energies between programs.

## Jobs

- Singlet `1247745_0` cancelled at the user's request after 58 min 28 s.
  Last accepted step: 48; energy: -2244.7410430222 Eh.
- Triplet `1247745_1` had already failed after 28 accepted steps with
  `ciah_macro_energy_rejected`; energy: -2244.790784197357 Eh.
- Spin/gradient diagnostic `1247926`: both spin checks passed. Endpoint
  conversion failed before the gradient test because the Molden shell header
  omits the optional shell scale. This is a diagnostic parser error, not an
  optimizer error. The parser was corrected; diagnostic `1247942` completed.

## Spin-space evidence

The current generated FULL basis has 32670 singlet structures and 50820 triplet
structures. These equal the spin-pure CAS(8,11) dimensions. Every generated
basis structure passes the explicit fermionic spin-raising test with zero
residual. This excludes a wrong-spin FULL structure generator for these inputs.
The reported target spin-square is not an independently measured expectation
value of the transformed CI carrier; carrier/action spin preservation remains
a separate check.

The orbital chart reports 37879 physical directions, equal to the core–active,
core–virtual, and active–virtual count for 89 core, 11 active, and 369 virtual
orbitals. Dimension agreement alone does not certify every direction or its
derivative.

## Confirmed termination mechanism

`optimization/backends/ciah.cpp` accumulates frozen accepted-point Hessian
images over microsteps. Real keyframes replace the working gradient but do
not replace that accumulated quadratic model. Macro acceptance then tests the
relaxed trial energy against the initial gradient and frozen accumulated
Hessian image using `optimization/neo/globalization.cpp`.

On rejection, it repeatedly halves the same direction. If the orbital step
becomes smaller than the existing minimum step, it returns
`ciah_macro_energy_rejected` and the outer optimizer exits. It does not refresh
the direction after this failure. The error name does not distinguish actual
energy increase, nonpositive predicted reduction, or inadequate actual/predicted
ratio. The original job did not record rejected-trial scalar diagnostics, so
none of those three numerical causes has yet been established for macro 29.

In particular, a keyframe-corrected accumulated step need not be a descent
direction for the initial gradient. Halving such a direction cannot guarantee
an acceptable sufficiently short step. This is a concrete model-consistency
risk, not evidence that changing an acceptance constant would fix the cause.

The current coefficient retraction is additive in raw orbital coefficients;
it is not a spherical exponential or an orbital rotation. The frozen chart's
linear gradient projection is appropriate to this particular retraction.
The finite CI normalization factor is present in the keyframe gradient.

## Energy discrepancy remains unresolved

Existing PySCF final energies are -2244.734377710447 Eh (singlet) and
-2244.782671264980 Eh (triplet). CPP's last accepted energies are lower by
6.665312 and 8.112932 mEh, respectively. Neither CPP run has converged.

The singlet trajectory passes close to the PySCF energy at step 23 and then
continues lowering the energy with a rising gradient. This demonstrates a
trajectory divergence, not a proof of either a correct lower minimum or an
incorrect derivative. The spin-generator certificate and the macro rejection
mechanism must not be used as an explanation for the energy difference without
additional derivative/carrier evidence.

Diagnostic `1247942` checks the analytic derivative along the negative reduced
gradient, normalized to unit physical norm, at the exported triplet endpoint
against two centered differences,
using the same current code and precision. The Molden export reconstructs an
approximate endpoint, not a bit-exact checkpoint of the failed macrostep.
The cancelled singlet did not export an endpoint orbital file.

## Completed directional derivative test

| Quantity | Value |
|---|---:|
| Analytic directional derivative | -0.0058411511060520883 |
| Centered difference, h = 0.001 | -0.0058392874962009955 |
| Absolute error, h = 0.001 | 1.863609851092754e-6 |
| Centered difference, h = 0.0003 | -0.0058409614211996095 |
| Absolute error, h = 0.0003 | 1.8968485247881223e-7 |
| Raw/chart derivative mismatch | 6.0715321659188248e-18 |

The derivative error decreases by about 9.8 times when h decreases by 3.33
times, consistent with centered-difference truncation convergence. This is
positive evidence for the relaxed gradient and its chart adjoint **in this
one direction**. It is not a certification of every orbital direction, the
coupled HVP, or spin preservation by the transformed carrier.

No optimizer thresholds or production source were modified. Next checks must
distinguish carrier spin preservation / physical variational space from
keyframe/coupled-model inconsistency, and record rejected macrostep scalars;
the current evidence does not establish the cause of the lower energies.

## Follow-up: basis, actual carrier, coupled derivative, and rejected step

Jobs `1247963` and `1247977` ran on Hanhai22. The sole production-source
change adds scalar `ciah_macro_trial` logging; it does not change numerical
decisions. A diagnostic-only `coupled-keyframe` mode compares the coupled
operator with the derivative of the same finite-normalized keyframe covector
used by production CIAH. It does not solve the relaxed structure response.

### Actual reference basis

The original PySCF checkpoint was copied from its existing Hanhai25 run.
For Fe, C, N, and H, its parsed shell exponents and contraction coefficients
equal `parse_gaussian.load` of the actual CPP input's `cc-pvdz.gbs`, with
`optimize=False`. The basis SHA256 is
`d53e7e02d55fb99cdf823d50cf13fec375e0cbe1220ff61816aee7ab38e38ec9`,
matching the original guess audit and repository file. The checkpoint has
37 atoms, 469 Cartesian AO functions, 186 electrons, zero charge, no spatial
symmetry, empty ECP, and no pseudopotential. Thus a differently named basis or
an ECP is not an explanation for this comparison.

### Actual transformed triplet carrier

Exporting `StructureAction::orthogonalize_structure_block` yields a
462-by-165 coefficient matrix. Its occupations are explicitly reordered into
PySCF's canonical string order before spin analysis; no assumed string order
or printed target spin-square is used.

- Carrier norm: 0.99999999999999434.
- PySCF measured spin-square: 1.999999999999978.
- PySCF spin-eigenvector residual norm: 2.487646800588102e-16.

There is no observed spin leakage in this actual triplet wavefunction. This
does not by itself certify every unprojected sigma image or every singlet
endpoint, which was not saved by the cancelled job.

### Coupled HVP versus production keyframe covector

The joint direction combines a physically normalized negative orbital
gradient with an independently chosen, structure-metric-normalized CI tangent.
The coefficient normalization derivative is included exactly as in production.

| h | Orbital image norm | Orbital error norm | Structure image norm | Structure error norm |
|---:|---:|---:|---:|---:|
| 0.001 | 46.6821156556 | 2.26217027e-5 | 11.5930118551 | 1.37279492e-5 |
| 0.0003 | 46.6821156556 | 2.03648550e-6 | 11.5930118551 | 1.23550313e-6 |

Both errors decrease by approximately 11.1 times, as expected for a centered
difference with the squared step-size ratio. The smaller-step relative vector
errors are approximately 4.36e-8 and 1.07e-7. This direction shows no coupled
HVP/keyframe inconsistency at the accepted point. It does not prove uniform
finite-step model accuracy throughout an entire macrostep.

### Bounded restart: concrete macro rejection

This restart uses the exported triplet endpoint, no retained seed, and a
two-macrostep cap. It is **not** an exact replay of the original macro 29.
The original orbitals/seed at the failed macro were not checkpointed.

At microstep 9, the frozen gradient prediction has norm 0.01480096, whereas
the actual keyframe gradient has norm 0.04762050. The mismatch is 0.04208605.
The current existing rule accepts that keyframe. Subsequent microsteps
produce an accumulated direction that is uphill for the initial gradient:

| Scale | Linear term | Quadratic term | Predicted reduction | Actual reduction |
|---:|---:|---:|---:|---:|
| 1 | 6.78865399e-5 | 1.40617575e-3 | -7.70974414e-4 | -1.95780881e-3 |
| 0.5 | 3.39432699e-5 | 3.51543937e-4 | -2.09715239e-4 | -2.89956079e-4 |
| 0.125 | 8.48581749e-6 | 2.19714961e-5 | -1.94715655e-5 | -2.00086943e-5 |
| 0.0009765625 | 6.62954491e-8 | 1.34103370e-9 | -6.69659660e-8 | -6.70088411e-8 |

Positive linear and quadratic terms imply negative predicted reduction for
**every positive scale**, not merely these sampled scales. Therefore halving
cannot satisfy the current acceptance contract. The diagnostic was stopped
after capturing this proof and repeated actual energy increases, avoiding
further pointless CI solves for ever-smaller scales. The completed coupled
derivative test took 99 s; the bounded restart was stopped after 264 s.

This establishes a failure of the accumulated microstep/globalization
combination at this restart point, not an explanation for previously accepted
lower energies. Changing the acceptance ratio cannot turn this direction into
a descent direction. The next algorithm investigation must address keyframe
anchoring, preservation of useful earlier candidates, and reconstruction of a
descent macrostep, not threshold tuning. Energy lowering versus PySCF remains
unresolved; it must not be attributed to spin contamination or ECPs after
these tests, nor assumed to be a correct lower minimum without further evidence.
