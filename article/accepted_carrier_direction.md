# Orbital directions from the accepted exterior carrier

Let $\Phi=\chi R$ be the accepted nonorthogonal active orbitals expressed in
an orthonormal active basis, with upper-triangular nonsingular $R$. Write
$T_\sigma=\wedge^{n_\sigma}R$ and let $X$ be the raw spin-string coefficients
of any selected VB expansion. The accepted carrier and its directional
derivative are

$$
Y=T_\alpha X T_\beta^T,\qquad
\dot Y=\dot T_\alpha X T_\beta^T
+T_\alpha X\dot T_\beta^T+T_\alpha\dot X T_\beta^T.
$$

The exterior representation is a group homomorphism. With
$K=\dot R R^{-1}$, its differential obeys

$$
\dot T_\sigma T_\sigma^{-1}
=d\Gamma_\sigma(K),\qquad
d\Gamma_\sigma(K)=\sum_{pq}K_{pq}a^\dagger_{p\sigma}a_{q\sigma}.
$$

Consequently the same derivative is obtained without replaying the accepted
exterior transformation:

$$
\boxed{\dot Y=d\Gamma_\alpha(K)Y+Yd\Gamma_\beta(K)^T
+T_\alpha\dot X T_\beta^T.}
$$

Compute $K$ by a small triangular solve, not an explicitly stored inverse.
Both $R$ and $\dot R$ are upper triangular. Their generated product support
is an upper-generator ideal, so this infinitesimal action stays on the same
exact symbolic support. No closure assumption on the selected VB structures
is made: the derivative is generated on the spin carrier and only subsequently
projected onto the original selected structures.

## Equal spin-exchange sectors

When identical spin spaces have $Y^T=sY$, $s=\pm1$, put
$G=d\Gamma_\alpha(K)Y$. Then

$$
Yd\Gamma_\beta(K)^T=sG^T,\qquad
\dot Y=P_s\left(2G+T_\alpha\dot X T_\beta^T\right),\quad
P_s Z=\tfrac12(Z+sZ^T).
$$

Thus only alpha generator links need be traversed. The already bounded input
tiles hold $2G$ (plus the transformed CI direction, when present) and are
projected in place. Unequal spin spaces use both generators and do not apply
this exchange identity.

## Storage and physical scope

Accepted compact $Y/(H-E)Y$ panels are already retained only when their whole
state fits the existing shared budget. A single directional component on the
same support fits within that same budget, so its tiles can be completed before
the in-place exchange projection. No extra complete carrier, larger workspace,
or thread-local coefficient panel is needed. If the accepted carrier is not
retained, the existing streamed differentiated exterior transform remains the
exact generator of the same directional data.

The one-body panel kernel uses one bounded shared dense input panel, borrowed
immutable beta spectator columns, and caller-owned output. It never forms a
dense spin-string generator matrix. Its work is linear in the number of acted-on
products times the single-replacement degree, rather than replaying an ordered
primal-and-tangent exterior jet. This is a reuse identity, not a reduction of
the product-space or Hamiltonian channel-GEMM asymptotic dimension.

Only generation of $\dot Y$ changes. The physical integral direction,
$H\dot Y+\dot H Y$, selected projection, representation derivative, integral
adjoint, and orbital HVP pullback are unchanged. In particular this is neither
a replacement of the nonorthogonal orbital update by a rotation nor removal of
the orbital-structure response.
