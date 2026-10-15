# Source-column contraction of arbitrary selected VB structures

The complete derivation, including the selected-space restriction proof,
Hamiltonian channel kernel, and implementation correspondence, is in
[the TeX note](source_column_hamiltonian.tex)
([compiled PDF](source_column_hamiltonian.pdf)).

## Exact expansion

Let $W$ map a selected VB coefficient vector $c$ to the original nonorthogonal
spin-product coefficients $X=\operatorname{reshape}(Wc)$. The selected
structures need not form a complete or transformation-invariant space.
Group the original nonzero support by its touched beta columns $J$:

$$
X=\sum_{j\in J} A_j e_j^T=A B^T,\qquad r=|J|.
$$

Here $A_j=X[:,j]$, and column $j$ of $B$ is the corresponding coordinate
vector. This is an exact source factorization, not an SVD, approximate rank
compression, or assertion that one VB structure has rank one. For example,
a product of six independent singlet pairs generally has 64 touched source
columns even when only one structure is selected.

For the exterior orbital maps $T_\alpha,T_\beta$, set

$$
L=T_\alpha A,\qquad R=T_\beta B,\qquad Y=LR^T.
$$

All subsequent contractions can be performed on $L$ and $R$, without
materializing the transformed product $Y$.

The current source plan improves the initial beta-column grouping by a minimum
mixed row/column cover of the original input support. In its bipartite graph,
each product is an edge between an alpha and a beta string. Assign every edge
to one covered endpoint, once only, giving the exact factorization

$$
X=\sum_{a\in C_\alpha}e_a u_a^T+
  \sum_{b\in C_\beta}v_b e_b^T=AB^T.
$$

Maximum matching certifies the minimum number of these coordinate-anchored
factors. Thus all formulas below also apply with
$r=|C_\alpha|+|C_\beta|\leq\min(|A_\Omega|,|J|)$.
This is not a numerical rank approximation. With a parity quotient the cover
uses representative input edges, while final sampling retains all original
output products. Section 3.2 of the TeX note gives the proof and an example
where a pure column grouping needs many sources but a mixed cover needs two.

If every selected structure has the same exact exchange parity $\eta$,
an additional source reduction is possible. Let $\Pi X=X^T$ and
$P_\eta=(I+\eta\Pi)/2$. Choose original product representatives
$a\leq b$ (strictly $a<b$ for odd parity), multiply off-diagonal
coefficients by two, and call the resulting matrix $X_{\mathrm{rep}}$.
Then $P_\eta X_{\mathrm{rep}}=X$. The transformed Hamiltonian and metric
preserve this sector, hence

$$
W^T T^T H T\operatorname{vec}(X_{\mathrm{rep}})
=W^T T^T H T\operatorname{vec}(X),\qquad T=T_\beta\otimes T_\alpha.
$$

The same holds with $H=I$ for the metric. Source factorization can therefore
use $X_{\mathrm{rep}}$, with only its touched source columns. The output
still samples **both** original exchange halves: the unprojected image of
$X_{\mathrm{rep}}$ is not a parity eigenvector. Even diagonal entries retain
weight one. Without a certified common parity the unreduced $X$ is used.

For a single covalent structure with six or seven singlet pairs, the source
counts are consequently 32 or 64, respectively, instead of 64 or 128.

## Hamiltonian factorization

Use the existing symmetric packed one-body channels
$E_P=E_{\alpha P}+E_{\beta P}$ and the symmetric kernel $K$ in

$$
H=\sum_{PQ}K_{PQ} E_P E_Q.
$$

The one-electron Hamiltonian and the contraction correction are already
absorbed into $K$. Consequently the following same-spin operators are
defined from this $K$; they should not separately be interpreted as physical
alpha/beta Hamiltonians with independently assigned one-electron energies.
Define

$$
D_\sigma=\sum_{PQ}K_{PQ} E_{\sigma P}E_{\sigma Q}.
$$

Commutativity of the two spin sectors and symmetry of $K$ give

$$
HY=(D_\alpha L)R^T+L(D_\beta R)^T
 +2\sum_{PQ}K_{PQ}(E_{\alpha P}L)(E_{\beta Q}R)^T.
$$

Both equal and unequal spin-string spaces obey this identity; no exchange
parity or selected-space closure is required.

Let

$$
\begin{aligned}
L_0&=T_\alpha^T L,&R_0&=T_\beta^T R,\\
L_H&=T_\alpha^T D_\alpha L,&R_H&=T_\beta^T D_\beta R,\\
L_P&=T_\alpha^T E_{\alpha P}L,&
V_P&=T_\beta^T\sum_Q K_{PQ}E_{\beta Q}R.
\end{aligned}
$$

Only the original touched products $(a,b)\in\Omega$ are sampled:

$$
\begin{aligned}
m_{ab}&=L_0[a,:]\cdot R_0[b,:],\\
h_{ab}&=L_H[a,:]\cdot R_0[b,:]
       +L_0[a,:]\cdot R_H[b,:]
       +2\sum_P L_P[a,:]\cdot V_P[b,:].
\end{aligned}
$$

The final images are $W^T\operatorname{vec}(h)$ and
$W^T\operatorname{vec}(m)$. They are exactly the restricted VB Hamiltonian
and metric actions. No assumption that the transformed vector remains in
the selected space enters the derivation.

## Causal exterior projection

The final sampling only needs rows in the original touched alpha and beta
sets. Computing all rows of $T_\sigma^T Z$ discards useful sparsity.
For an elementary exterior shear with coefficient $s$ and sign $\epsilon$,
the adjoint update is

$$
z_i\leftarrow z_i+s\epsilon z_j.
$$

If row $i$ is not needed by any subsequent update or final sample, this edge
has no effect on the answer. Starting from the requested final rows, traverse
the adjoint operations in reverse chronological order: retain an edge exactly
when its updated row is needed, and then mark its source row as needed.
The retained graph is a causal dependency closure, independent of coefficient
magnitudes. Execute the retained edges in their original adjoint order.

Writing $P_A$ for selection of the required rows, the operation remains
$P_A T_\sigma^T Z$ exactly. It is not a projection of $Z$ onto the original
selected VB subspace. Intermediate rows outside $A$ are retained whenever a
path from them reaches $A$. Consequently arbitrary selected spaces require
no closure assumption here either.

The plan stores shear-edge indices and one spin-row bitmap, not a compound
matrix or product-space table. It is generated once per source contraction
and shared by all source panels, channels and RHS vectors.

## Source and channel tiles

Split the source columns into panels of width $s$ and channels into panels
of width $p$. For each source panel and channel panel $P$, construct
$E_{\alpha P}L$ and $E_{\beta P}R$. Accumulate the weighted channel panels
by streaming $Q$ and multiplying by $K_{PQ}$. Gather their same-spin
contributions immediately, pull back the cross-spin panels and sample
$\Omega$, then release the panel. All source panels contribute to the same
small original-support output arrays.

The workspace is

$$
O\big((U_\alpha+U_\beta)s p\big),
$$

with no additional $U_\alpha U_\beta$ array and no per-thread carrier.
The implementation uses the same 64 MiB numerical-workspace ceiling as the
product-tile sigma. A smaller explicitly requested source workspace reduces
$s$ and $p$; it never changes the represented VB space.

Channel and source slabs are allocated once for an action and reused across
source panels, channel panels and all RHS. Full-row views represent partial
panels without resizing their backing allocations. This is required to keep
streaming costs proportional to reading and updating the slabs, rather than
repeatedly mapping and first-touching their pages.

## Work and choice of contraction order

Let $P=n(n+1)/2$ be the number of channels, $t=|\Omega|$, and $Q_\sigma$
the number of single-spin exterior shear pairs. With a full channel panel,
the leading source-order work per RHS is

$$
O\left[r(U_\alpha+U_\beta)P^2
 +rP(Q_\alpha+Q_\beta)+trP\right].
$$

The generated product-tile order instead has channel-GEMM work
$O(N_{\mathrm{closure}}P^2)$ plus exterior transforms and projection. If
$r\ll\min(U_\alpha,U_\beta)$ and the transformed closure is broad, the
source order removes a factor comparable to $U/r$ from the leading carrier
work. When many original source columns are touched, this advantage can
disappear. In particular, a complete covalent selected space can have $r=U$.
No universal linear dependence on the number of selected structures follows.

The implementation compares arithmetic estimates for these two exact
contraction orders. The source estimate includes the two channel GEMMs,
spin links, shear pullbacks and support sampling. The product order is first
bounded by the complete carrier; if necessary, its generated closure is
counted only until its channel-GEMM work exceeds the source estimate. This
planning happens once per immutable action, not once per Davidson vector.
It contains no fitted threshold, memory-based physical admission, or
approximate rank selection. Arithmetic estimates need native measurements
to establish a wall-time benefit.

Channel tiling can repeat spin-link scatters between channel panels, but it
does not alter the two-body formula. Only source columns and channel blocks
are stored; discarded panels are consumed immediately.

This implementation targets the fixed-point H/M action. Integral adjoints,
orbital gradients and HVP are still evaluated by their existing exact kernels;
the source-order scaling is not claimed for those consumers yet.
