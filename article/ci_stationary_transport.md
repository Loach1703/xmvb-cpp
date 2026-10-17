# CI stationarity and a residual-balanced transported CIAH model

## Physical problem and finite transport

Let $p$ denote the nonredundant orbital coordinates measured from the accepted
point, and let $q$ denote the independent structure coordinates at that point.
The fixed linear map $Q$ embeds $q$ into structure coefficients, including the
existing state-weight scaling and horizontal projection. Define

$$
A(q)=C_0+Qq,\qquad F(p)=\bigl(T_\alpha(p)\otimes T_\beta(p)\bigr)W,
\qquad M(p)=F(p)^TF(p).
$$

$W$ maps the **selected VB structures** into the original touched spin-string
products. Its columns need not span a closed space under orbital changes.
The tensor-product notation defines the map, not a stored product matrix.
The production implementation applies it by bounded factor/tile contractions.

The current finite CI transport is the selected-space least-squares solution

$$
\tau(p,q)=M(p)^{-1}F(p)^T
\left[F(0)A(q)+DF(0)[p]C_0\right].
$$

For the positive selected metric, its first jet is

$$
\tau(0,0)=C_0,\qquad D_p\tau(0,0)=0,\qquad D_q\tau(0,0)=Q.
$$

This follows by differentiating the normal equation. The orbital derivative
of its right-hand side equals $DM(0)[p]C_0$, leaving $M(0)D_p\tau=0$.
No ambient-space closure or full-CAS identity is invoked.

Let $\mathcal E(p,C)$ be the physical, normalized single-state or equal-weight
state-average energy, including inactive and nuclear terms. The physical
transported energy and its affine-coordinate counterpart are

$$
E_\tau(p,q)=\mathcal E(p,\tau(p,q)),\qquad
E_A(p,q)=\mathcal E(p,A(q)).
$$

Their gradients agree at the origin. The chain rule gives, for arbitrary
directions $u$ and $v$ in the combined orbital/structure coordinates,

$$
D^2E_\tau(0)[u,v]=D^2E_A(0)[u,v]
+\left\langle g_{c,0},D^2\tau(0)[u,v]\right\rangle_F,
\qquad g_{c,0}=\nabla_C\mathcal E(0,C_0).
$$

The last term is a **coordinate curvature**, not the physical orbital--CI
response. With $C_0^TM_0C_0=I$ and common state weight $w$, the implemented
origin gradient is

$$
g_{c,0}=2w\left(H_0C_0-M_0C_0\Lambda_0\right),\qquad
\Lambda_0=C_0^TH_0C_0.
$$

At an exact selected-space CI solution $H_0C_0=M_0C_0\Lambda_0$, this gradient
and its coordinate-curvature contraction vanish. For an inexact numerical
eigenpair they need not vanish. Simply replacing them by zero is not an exact
algebraic elimination and is **not** the implementation proposed here.

## Residual-balanced model: exact cancellation also for inexact CI

Freeze the actual computed $g_{c,0}$ once per accepted macro point and define
the internal CIAH model

$$
\widehat E(p,q)=E_\tau(p,q)
-\left\langle g_{c,0},\tau(p,q)-A(q)\right\rangle_F.
$$

The coefficient of this counterterm is fixed by the chain rule; there is no
user parameter, residual threshold or conditional path. The first-jet
identity implies, **for any finite origin CI residual**,

$$
\widehat E(0)=E_\tau(0),\qquad
D\widehat E(0)=DE_\tau(0)=DE_A(0),\qquad
D^2\widehat E(0)=D^2E_A(0).
$$

Thus the numerical CI residual remains in the CIAH right-hand side. Only its
contraction with the nonlinear transport's second jet is cancelled. The
expensive transport-curvature operator and its origin adjoint preparation
are unnecessary for this model, including when $g_{c,0}$ is deliberately large.

This is an explicitly defined internal model, **not** a claim that the
physical transported Hessian equals the affine Hessian away from CI
stationarity. The physical integral, gradient and HVP kernels are unchanged.
At exact CI stationarity $g_{c,0}=0$, the model is the physical transported
energy at every finite $(p,q)$, not just through second order. As numerical CI
residuals vanish, the difference vanishes continuously without a cutoff.

## Consistent finite keyframes

The finite keyframe must use the energy and gradient of the **same** model.
Writing $g_C=\nabla_C\mathcal E(p,\tau)$ and $g_p$ for the fixed-$C$ orbital
gradient, differentiation yields

$$
\nabla_p\widehat E=g_p+(D_p\tau)^T(g_C-g_{c,0}),
$$

$$
\nabla_q\widehat E=(D_q\tau)^T(g_C-g_{c,0})+Q^Tg_{c,0}.
$$

The existing implicit transport pullback therefore uses the right-hand side
$g_C-g_{c,0}$, rather than $g_C$, in the **same one metric solve**. The final
structure gradient adds the frozen affine term $Q^Tg_{c,0}$. The model energy
subtracts the matching scalar counterterm. This neither introduces another
metric inverse nor retains any new spin-product array.

All finite nonclosure/residual terms in the existing transport pullback are
retained. Setting $D_p\tau$ to zero away from the origin would be incorrect.
Likewise, using a corrected Hessian with an uncorrected keyframe model would
break derivative consistency; that mixed implementation is excluded.

## Physical coupled response remains

At the origin the model action is the existing physical affine-coordinate
joint action,

$$
\begin{pmatrix}y_o\\y_c\end{pmatrix}=
\begin{pmatrix}
\mathcal A(p)+\mathcal B^\dagger(q)\\
\mathcal B(p)+\mathcal C(q)
\end{pmatrix}.
$$

These are matrix-free operator actions, not assembled Hessian blocks.
Orbital integral response, the orbital-to-CI forcing, CI Hessian/metric
actions, the CI-to-orbital integral adjoint and the coarse preconditioner all
remain. In particular, no relaxed physical response or excited-state
contribution is removed by eliminating the coordinate counterterm.

Macro acceptance and completion still evaluate $\mathcal E$ and solve the
physical CI problem. Outer energies, gradients, stopping criteria and
acceptance tests remain physical; only the internal CIAH keyframe
energy/gradient is $\widehat E$. Diagnostics distinguish physical and model
energies. No new convergence or stability test is added.

## Implementation and validation contract

1. `CiTransport` retains its finite selected-space retraction and one-solve
   pullback, exposes a residual-balanced model evaluation, and keeps the
   actual origin CI residual for the initial right-hand side.
2. Delete the separate residual-weighted origin curvature, its cached metric
   inverse image/adjoints, and the now-unused fused transport-curvature sweep.
   Generic cross-overlap actions, adjoints and finite transport gradients
   remain because they have physical consumers.
3. Test energy/gradient finite differences on arbitrary nonclosed selected
   spaces with **large nonstationary** origin residuals. Compare origin HVPs
   of the balanced model with affine-coordinate HVPs. Include a stationary
   generalized-eigenvector case and equal-weight multistate cases. The
   nonstationary fixture prevents an implementation based on a smallness gate
   from passing.
4. Preserve existing exact HVP/coupled-response regressions and run only the
   new C12/LOFLEA version on Hanhai22, 32 cores, to convergence. Compare against
   already saved energies, steps, HVP counts, native wall time and RSS.

No speedup or preserved numerical convergence is inferred from this derivation;
the production validation must establish them separately.
