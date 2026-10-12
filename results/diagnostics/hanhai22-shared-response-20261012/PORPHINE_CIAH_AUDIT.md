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
