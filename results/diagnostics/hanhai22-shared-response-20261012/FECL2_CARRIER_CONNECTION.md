# Carrier connection in the coupled nonorthogonal model

## Status and purpose

This is a derived diagnostic hypothesis, not an established implementation
bug or an implemented optimizer change. The FeCl2 tail residual is mostly
orbital; neither orbital-block Jacobi nor metric-orthogonal AH vectors
reduced the work. The next test should separate physical integral coupling
from coupling caused by changing the representation of the CI coefficients.
No orbital direction, spin space, retraction or threshold should be changed
before this decomposition is measured and checked.

## Structure and spin-string spaces

At orbital coefficients $b$, let $F(b)$ map VB structure coefficients
$\mathbf C$ into the globally orthogonal spin-string product carrier:

$$
y=F(b)\mathbf C,\qquad M=F^TF,\qquad H=F^T h F.
$$

Here $h$ is the orthogonal-carrier Hamiltonian, $M$ and $H$ are structure
overlap and Hamiltonian matrices, and $y$ is a flattened spin-string product
tensor. These formulas do not introduce a determinant intermediate in the
implementation. Neither $F$ nor a complete dense $H$ or $M$ is proposed as
a production allocation.

For FeCl2, $U_\alpha=\binom87=8$ and $U_\beta=\binom83=56$ give 448
carrier entries. The complete target-spin structure space has dimension
$448-\binom88\binom82=420$, not 448; the extra carrier spin sector must not
be mistaken for an additional variational freedom.

## Exact mixed derivative

Write $r=(h-E)y$ at $\mathbf C^TM\mathbf C=1$. For an orbital direction
$p$, denote $\delta F=DF[p]$ and $\delta h=Dh[p]$. Before projecting onto
the structure tangent, the mixed derivative of the coefficient gradient is

$$
Bp=2\delta F^Tr+
2F^T(h-E)\delta F\mathbf C+
2F^T\delta h\,y-
2\delta E F^Ty.
$$

The last term vanishes after projection onto an overlap-horizontal
structure tangent, because its columns $Z$ satisfy $Z^TM\mathbf C=0$.
In a complete invariant target-spin space, $r=0$ at the relaxed eigenstate.
Then only two terms remain:

$$
Z^TBp=2Z^TF^T(h-E)\delta F\mathbf C+
2Z^TF^T\delta h\,y.
$$

The first is a representation (carrier connection) term, and the second
contains the physical integral response. For an arbitrary selected VB
subspace, $F^Tr=0$ does **not** imply $r=0$. The $\delta F^Tr$ term and any
closure defect must therefore be retained. Completeness must not be inferred
merely from complete unique-spin-string lists.

## Connection direction and its closure certificate

Remove normalization from the basis-induced carrier displacement:

$$
d_y=(I-yy^T)\delta F\mathbf C.
$$

For the unscaled tangent map $F_Z=FZ$, find its metric least-squares
coordinates from

$$
(F_Z^TF_Z)\gamma=F_Z^Td_y,\qquad
e_{\mathrm{closure}}=d_y-F_Z\gamma.
$$

This is an overlap-metric solve, **not** a structure-Hessian response solve.
When the selected structure space is invariant under the spin-scalar
one-body generator and complete in the target spin, the closure defect is
zero in exact arithmetic. Otherwise only a projected connection is obtained;
the remaining carrier direction is physical and must not be discarded.
In code, tangent coordinates are scaled by $\sqrt{2w}$; this factor must be
included explicitly when relating $\gamma$ to `StructureTangentOperator`.

If $K=2F_Z^T(h-E)F_Z$ is the unscaled structure curvature, completeness gives

$$
Z^TBp=K\gamma+2Z^TF^T\delta h\,y.
$$

Thus the coupled forcing can be large even when the physical integral
coupling is modest. This is a possibility to measure, not the conclusion
of the current tests.

## A possible preconditioner, not a new variational space

If the direction-dependent connection defines $\gamma=\Gamma p$, use the
linear shear

$$
\begin{pmatrix}p\\q\end{pmatrix}
=J\begin{pmatrix}p\\s\end{pmatrix},\qquad
J=\begin{pmatrix}I&0\\-\Gamma&I\end{pmatrix}.
$$

The transformed Hessian is a congruence, not the deletion of a coupling:

$$
J^T\begin{pmatrix}A&B^T\\B&K\end{pmatrix}J
=\begin{pmatrix}
A-B^T\Gamma-\Gamma^TB+\Gamma^TK\Gamma & B^T-\Gamma^TK\\
B-K\Gamma & K
\end{pmatrix}.
$$

The full AH metric must transform as $J^TGJ$ as well. Alternatively this
shear may be used only to construct a preconditioner, leaving the real AH
operator and metric unchanged. A block-diagonal metric in the new coordinates
cannot be substituted without acknowledging that it changes the AH model.
No performance benefit is assumed: applying $\Gamma$ and its adjoint may
cost enough to outweigh better conditioning.

## Existing code and next bounded diagnostic

`StructureAction::apply_integral_direction` already forms the relative
upper-triangular orbital-transform direction, then
`sigma.apply_one_body_generator(spin_vectors, relative_orbital_direction)`.
This is the basis-induced $\delta F\mathbf C$, distinct from a physical CI
response. Its sigma action and the integral-direction sigma action are
currently added before projection back to structures. The decomposition can
therefore be observed at this point without an additional orbital HVP or
AO-integral traversal.

The diagnostic should report, in the complete structure dual norm:

1. Total coupled forcing, basis-induced forcing and physical integral forcing.
2. Their cross pairing; separate norms alone can hide strong cancellation.
3. The metric solve's unchanged true residual and the carrier closure defect.
4. Reconstruction of total forcing from these terms, with the current
   normalization, state weight and tangent-coordinate factors.

Use FeCl2 first and existing accepted-point data. Do not infer a bug from
large components alone, form a dense coupling matrix, extend the selected
space, or start another iron-porphine long run before this check.
