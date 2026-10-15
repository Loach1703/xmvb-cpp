# Exact source-factor contraction for selected-space CI transport

## Representation and scope

Let $W$ be the fixed map from selected VB structure coefficients to the original
alpha/beta spin-string products, and define $X(c)=\operatorname{reshape}(Wc)$.
The orbital chart $R$ is upper triangular, with the exterior maps

$$
T_\sigma(R)=\wedge^{n_\sigma}R,\qquad
F(R)c=T_\alpha(R)X(c)T_\beta(R)^T.
$$

The selected structures can be arbitrary and need not be closed under orbital
transformations. Alpha and beta string spaces are treated independently.
This derivation concerns cross overlaps, CI transport, and their orbital
pullbacks, not the complete Hamiltonian sigma or integral adjoint.

## Exact source factorization

Let $\mathcal J$ contain the beta columns touched by the original selected
structure expansion. Partition it into panels $J$ and let $E_J$ consist of the
corresponding columns of the beta identity operator. Then, exactly,

$$
X=\sum_J X_JE_J^T,\qquad X_J=XE_J.
$$

This identity uses structural zeros only. There is no singular-value truncation,
numerical rank threshold, or approximation. Its source width
$r=|\mathcal J|$ is not necessarily the minimal algebraic rank of $X$.

For two orbital charts, define the single-spin cross operators

$$
K_\sigma=T_{\sigma,L}^TT_{\sigma,R}.
$$

Associativity gives

$$
F_L^TF_Rc
=W^T\operatorname{vec}\left[
\sum_J (K_\alpha X_J)(K_\beta E_J)^T
\right].
$$

No $K_\sigma$ matrix is formed: its factors are applied successively.
For each panel, set $A_J=K_\alpha X_J$ and $B_J=K_\beta E_J$.
At an original touched product $(a,b)$, accumulate

$$
z_{ab}\mathrel{+}=A_J[a,:]\cdot B_J[b,:].
$$

Finally apply $W^T$ once. Crucially, neither the orthogonal product tensor nor
its reachable product support is generated. The sum is evaluated only at the
original touched products required by $W^T$.

This replaces the earlier output-panel implementation, which still generated
panels over the transformed support and compound beta rows. That implementation
has been removed from the transport path rather than retained as a fallback.

## Directional action

At fixed structure coefficients, differentiate the same factors:

$$
\dot A_J=\dot K_\alpha X_J,\qquad
\dot B_J=\dot K_\beta E_J,\qquad
\dot K_\sigma
=\dot T_{\sigma,L}^TT_{\sigma,R}
+T_{\sigma,L}^T\dot T_{\sigma,R}.
$$

The sampled output tangent is

$$
\dot z_{ab}\mathrel{+}=
\dot A_J[a,:]\cdot B_J[b,:]
+A_J[a,:]\cdot\dot B_J[b,:].
$$

Thus the action and its derivative have the same source width and panel
storage order. Neither needs a generator halo in the product space.

## Reverse contraction

For a scalar $f=c_L^TF_L^TF_Rc_R$, write $X_L=X(c_L)$ and use the source
panels of $X_R=X(c_R)$. Define

$$
U_J=T_{\alpha,R}X_{R,J},\qquad
V_J=T_{\beta,R}E_J,\qquad
A_J=T_{\alpha,L}^TU_J,\qquad
B_J=T_{\beta,L}^TV_J.
$$

The seeds for the two cross factors are

$$
P_{\alpha,J}=X_LB_J,\qquad
P_{\beta,J}=X_L^TA_J.
$$

Both are sparse matrix–dense panel products over the original touched
spin-string products. Form

$$
Z_{\alpha,J}=T_{\alpha,L}P_{\alpha,J},\qquad
Z_{\beta,J}=T_{\beta,L}P_{\beta,J}.
$$

For a spin-string sector, define the generator adjoint

$$
[G_\sigma(P,Q)]_{pq}
=\langle P,E_{\sigma,pq}Q\rangle_F,
$$

where $E_{\sigma,pq}$ is the one-body excitation generator acting on rows.
Since

$$
\dot T_\sigma=E_\sigma(\dot R R^{-1})T_\sigma,
$$

the two chart gradients are

$$
\nabla_{R_L}f=
\operatorname{upper}\left\{
\sum_J[
G_\alpha(U_J,Z_{\alpha,J})
+G_\beta(V_J,Z_{\beta,J})
]R_L^{-T}
\right\},
$$

$$
\nabla_{R_R}f=
\operatorname{upper}\left\{
\sum_J[
G_\alpha(Z_{\alpha,J},U_J)
+G_\beta(Z_{\beta,J},V_J)
]R_R^{-T}
\right\}.
$$

Only upper generator entries are needed for an upper-triangular chart.
The existing single-spin excitation topology supplies them; there is no
product-space excitation graph or transition-density array.

## Pullback derivative

All of $U,V,A,B,P,Z$ are differentiated, including both inputs of every
bilinear contraction:

$$
\dot G(P,Q)=G(\dot P,Q)+G(P,\dot Q).
$$

With $I=R^{-T}$,

$$
\dot I=-I\dot R^TI,\qquad
\dot{\nabla_R f}
=\operatorname{upper}(\dot G I+G\dot I).
$$

This retains the derivative of the inverse chart and all changes of the
seeds. It is the derivative of the same scalar cross overlap, not a
frozen-chart or frozen-response approximation. Structure-coefficient
directions are supplied separately by the CI-transport caller.

## Work and memory

Let $s$ be the number of original touched products, $r$ the number of touched
beta strings, $L_\sigma$ the number of single-spin elementary-transform links,
and $b$ the source-panel width. Per RHS, the action, its directional derivative,
and its pullbacks require, up to different constant factors,

$$
O\!\left[
r(L_\alpha+L_\beta+U_\alpha+U_\beta+s)
+\operatorname{nnz}(W)
\right].
$$

Using $L_\sigma=O(n_{\rm act}^2U_\sigma)$ gives

$$
O\!\left[
r\{n_{\rm act}^2(U_\alpha+U_\beta)+s\}
+\operatorname{nnz}(W)
\right].
$$

Panel workspace is $O(b(U_\alpha+U_\beta))$; touched coefficients and indexing
use $O(s+U_\alpha+U_\beta)$ in addition to the existing expansion and
single-spin topology. Threads share the factors. No thread receives its own
product tensor.

The byte budget bounds panel temporaries, not the caller-owned vectors and
immutable topology. It changes batching only, never the physical space.
This chain does not use the compound-row cache.

For one CAS(12,12) covalent structure, $r=64$ while $U_\alpha=U_\beta=924$.
The transform's source-column count is therefore reduced by a factor of
$924/64=14.4375$ relative to a full-column transform; this is an operation-count
comparison, not an asserted wall-time speedup. For all 132 covalent structures,
$r=924$, so this particular source-width advantage disappears. The elimination
of output-product traversal and compound-row generation still applies.

This is not a polynomial-scaling solution to arbitrary full active spaces.
For a full covalent family with complementary alpha/beta occupations,
$s=U_\beta$ and $r=U_\beta$, so the work remains quadratic in string counts
(apart from orbital factors). For dense original product support,
$s=U_\alpha U_\beta$, the sampled factor products can instead cost
$O(U_\alpha U_\beta^2)$. The restricted metric solve and separate
Hamiltonian/response kernels also retain their own costs.

## Validation

Small-system tests compare cross-Gram actions to independent occupied minors,
gradients to directional finite differences, and pullback directions to
finite differences of the gradient. Cases include arbitrary nonclosed
selected spaces, sparse supports, unequal spins, identity charts, and narrow
versus wide source panels.

CAS(12,12) fixed-orbital kernel tests use 1, 3, 20, and 132 covalent structures.
Production regression uses 240 COV, ION(0:1), ION(0:2), ION(0:3), and FULL with
unchanged inputs, optimizer, and tolerances. Timings and RSS are reported
separately from the mathematical operation count; historical binaries are
not rerun.
