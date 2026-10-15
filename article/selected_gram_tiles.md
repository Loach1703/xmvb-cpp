# Selected-space Gram actions and differentiated tile contractions

## Scope

This derivation concerns CI transport and its orbital pullbacks for an arbitrary
selected VB structure space. It does not replace the Hamiltonian sigma action,
change the selected structures, or assume closure of those structures under an
orbital transformation. All coefficients and orbital charts are real.

Let $W$ map structure coefficients to touched alpha/beta spin-string products.
For an upper-triangular orbital chart $R$, write

$$
F(R)c = T_\alpha(R)X(c)T_\beta(R)^T,
\qquad X(c)=\operatorname{reshape}(Wc),
\qquad T_\sigma(R)=\wedge^{n_\sigma}R.
$$

The notation $X$ specifies a linear map, not an instruction to allocate a full
product array. In the implementation, its coefficients are accumulated only at
the touched product indices. Alpha and beta need not have the same electron
number or the same string lists.

## Exact panel factorization

For a panel $J$ of output beta strings, define

$$
B_{\sigma,J}=T_\beta(R_\sigma)[J,:],\qquad
Y_{\sigma,J}=T_\alpha(R_\sigma)X_\sigma B_{\sigma,J}^T,
\qquad \sigma\in\{L,R\}.
$$

Here $L,R$ label the two charts, not the electron spins. Then

$$
c_L^T F(R_L)^T F(R_R)c_R
=\sum_J \langle Y_{L,J},Y_{R,J}\rangle_F.
$$

For an action rather than a scalar contraction, each right image is pulled
through the left alpha transform and sampled at the original touched products.
All panel contributions are accumulated before applying $W^T$ once.
No selected-space invariance has been used.

Upper-triangular transforms only move occupied orbital indices downward. Their
directional derivatives have the same structural reachability, including at
charts whose off-diagonal numerical entries vanish. The support planner must
therefore use structural reachability, not a numerical magnitude cutoff.

## Directional action

The panel tangent is

$$
\dot Y_J=\dot T_\alpha X B_J^T+T_\alpha X\dot B_J^T.
$$

Consequently,

$$
\frac{d}{dt}(F_L^TF_R)c
=\dot F_L^TF_Rc+F_L^T\dot F_Rc.
$$

The primal and tangent use the same panel sweep and exact support. A primal
structure projection is not computed when only the tangent output is needed.

## Pullbacks without a full product carrier

Let $E_\sigma(K)$ denote the one-body generator on the spin-string sector.
The exterior representation obeys both differential identities

$$
\dot T_\sigma=E_\sigma(\dot R R^{-1})T_\sigma
=T_\sigma E_\sigma(R^{-1}\dot R).
$$

For an alpha panel, the first identity gives a generator contraction over its
rows. If $G_\alpha$ is its adjoint with respect to $\dot R R^{-1}$, the chart
gradient is $G_\alpha R^{-T}$.

For the beta rows, first form the seed using the sparse original coefficients:

$$
Q_{L,J}=\left[T_\alpha(R_L)^T Y_{R,J}\right]^T X_L.
$$

The second identity yields

$$
\langle Q_{L,J},\dot B_{L,J}\rangle_F
=\langle Q_{L,J},B_{L,J}E_\beta(R_L^{-1}\dot R_L)\rangle_F.
$$

If $G_\beta$ is the generator adjoint in this expression, its chart gradient
is $R_L^{-T}G_\beta$. Sum the alpha and beta contributions and retain only
upper-triangular entries. Interchanging left and right gives the other chart
gradient. This does not require neighboring beta output panels or a complete
product-space transition density.

Only upper one-body generators contribute: in the upper part of
$G_\alpha R^{-T}$ and $R^{-T}G_\beta$, every accessed generator entry is upper
triangular. The existing single-spin excitation topology supplies these
generators; no product-space excitation table is built.

Differentiating these same contractions gives the pullback tangent. With
$I=R^{-T}$,

$$
\dot I=-I\dot R^T I,
\qquad
\dot g_\alpha=\dot G_\alpha I+G_\alpha\dot I,
\qquad
\dot g_\beta=\dot I G_\beta+I\dot G_\beta.
$$

Each $\dot G$ is obtained by differentiating both bilinear inputs. Thus the
HVP connection retains the derivative of the inverse chart as well as both
image tangents; it is not a frozen-chart approximation.

## Storage and scope of the speed claim

For $s$ touched original products and panel width $b$, the per-RHS action uses
$O(s)$ coefficients and outputs and $O(b(U_\alpha+U_\beta))$ panel workspace,
in addition to the structure vectors, per-spin excitation topology, and
bounded compound-row caches. The pullback has the same panel storage order.
There is no mandatory full $U_\alpha U_\beta$ coefficient image in this chain.
For full original support, $s$ itself can equal $U_\alpha U_\beta$; this method
does not claim to remove the physical full-space scaling.

The workspace argument bounds panels, not all input/output storage or all
immutable caches. Caches retain the existing fixed byte limit and are shared
across threads. A thread may allocate a single-spin support mask, not its own
product array.

This is a support-aware factorization and derivative reformulation, not a proof
that total SCF time scales linearly with the number of selected structures.
Beta compound-row generation, generator contractions, AO operations, and
metric iteration counts can still dominate. Native timings must establish
the wall-time benefit separately.

## Validation

The unit test uses an independent occupied-minor oracle for cross-Gram actions,
finite differences and adjoint duality for first and second derivatives,
unequal spin sectors, arbitrary nonclosed selected spaces, identity charts
with vanishing minors, sparse selected supports, and single-column versus
wide-panel execution. Production validation uses the existing 240 COV,
ION(0:1), ION(0:2), ION(0:3), and FULL inputs without changing convergence
tolerances or rerunning historical executables.
