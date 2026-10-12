# Natural-active-frame analytic diagonal: controlled preconditioner experiment

Rejected experiment; none of the natural-frame implementation is retained.

## Motivation

The accepted diagonal currently uses a Cholesky-orthonormalized active
frame. This is a valid but otherwise arbitrary orientation. A diagonal
approximation is not invariant under active basis mixing, even though
the full variational equations are. Active 1-RDM off-diagonal components
therefore remain outside this particular diagonal approximation.
Test a representation chosen by the accepted density, not a new threshold
or molecule-specific parameter.

## Coordinate transformations

Let Phi denote the unchanged physical active orbital matrix and
R^T R=Phi^T S Phi. The existing orthogonal frame is Q=Phi R^{-1}.
For accepted contravariant active 1-RDM gamma, form

$$
\gamma_Q=R\gamma R^T,
\qquad \gamma_Q=O n O^T,
\qquad Q_{\mathrm{nat}}=QO.
$$

No occupation is truncated. Natural occupations may include zero values;
no occupation floor or active-space selection is introduced. The new
source transforms are

$$
R_{\mathrm{nat}}=O^TR,\qquad
R_{\mathrm{nat}}^{-1}=R^{-1}O.
$$

The existing 1-/2-RDM and active-pair integral transformations use these
matrices consistently. Thus the analytic diagonal is recalculated in
the natural frame rather than merely rotating the old diagonal (which
would not be the diagonal of the transformed exact operator).
The physical orbitals, structure coefficients and spin-string carrier
are not rotated. No NEO coordinate, gradient, HVP or physical metric
changes. The full-AO inverse still uses its two contravariant source maps;
HAO still pulls the same analytic model into its strictly local tangent.

## Preserve the existing step limit

The component cap currently depends on the analytic frame, so simply
rotating that frame would inadvertently change globalization too.
The experiment explicitly restores components to the original frame
before evaluating the maximum. For full OEO, with core-to-active rows
and external-to-active columns respectively,

$$
\kappa_{ac}^{\mathrm{old}}=O\kappa_{ac}^{\mathrm{nat}},
\qquad
\kappa_{va}^{\mathrm{old}}=\kappa_{va}^{\mathrm{nat}}O^T.
$$

HAO restores both the active rows of each local external-component vector
and the active source weights before applying its existing component
bound. Only the n_active-by-n_active rotation is stored, not a second
AO frame or a per-thread array. This isolates preconditioner orientation
from changes to step limits, supports or retraction.

## Validation and rejected numerical result

Independent test references transform full four-index AO integrals and
the original density tensors to the reported natural frame. Every analytic
diagonal entry is checked for exact/RI, signed and zero/dependent sources,
and zero/one/two core orbitals. Separate HAO/OEO chart tests verify the
same nonredundant dimension, identical physical update and original-frame
maximum component under an explicit nontrivial active rotation.

Hanhai22 1249642 completed on anode31. All five unit executables passed,
including the independent natural-frame diagonal and original-frame cap
checks, followed by normal original-guess FeCl2 convergence under the
unchanged existing criteria. Nevertheless its endpoint is unacceptable:

| Variant / job | Macros | HVPs | Energy / Eh | GNU wall / s | GNU RSS / MiB |
|---|---:|---:|---:|---:|---:|
| Retained frame, archived 1249596 | 12 | 251 | -2181.624464871899 | 8.30 | 516.23 |
| Natural frame, 1249642 | 11 | 142 | -2181.624176749723 | 6.25 | 503.05 |

The natural-frame endpoint is +0.288134010 mEh above the independent
PySCF coupled-Newton reference. Its last reported physical orbital
gradient norm is 8.28045e-4 and adjacent accepted energy change is
7.01348e-8 Eh. Passing the shared stopping criteria does not prove that
this endpoint meets the requested reference-energy accuracy. Do not
claim that the action reduction makes this a successful optimizer change.
Do not label the endpoint a distinct local minimum without a curvature
diagnostic: it may still lie on a very soft unresolved direction.

The unit evidence supports the designed preconditioner transformation,
unchanged nonredundant update and cap; it does not prove improved nonlinear
convergence. This is exactly why an actual molecular test is required.
The experimental natural-frame fields, transform, cap restoration and
altered unit fixtures are removed. The original source and references
are restored; no selector, compatibility code or fallback remains.
No iron long run follows this rejected experiment.
