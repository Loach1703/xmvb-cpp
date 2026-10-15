# Exact spin-exchange quotient for projected direct-CI contraction

## Scope and representation

Let $W$ expand a chosen set of VB structures into alpha/beta spin-string
products. The selected set need not be complete or closed under orbital
transformations. Suppose the alpha and beta string lists coincide, and define
the exchange operator by

$$
\Pi X=X^T.
$$

The condition required here is stronger than equal spin-string dimensions:
every column of $W$ must have the same exchange parity $\eta=+1$ or $-1$,

$$
\Pi W=\eta W.
$$

This is a statement about the actual expansion coefficients, not an inference
from $M_S=0$, nor a total-spin projection. A mixed-parity space does not meet
the condition. The implementation compares the sparse coefficients at
$(a,b)$ and $(b,a)$ exactly, once when the structure map is built; it neither
rounds coefficients nor introduces a numerical symmetry threshold.

For identical spatial orbital transformations in both spin sectors,

$$
F c=T X(c)T^T,\qquad X(c)=\operatorname{reshape}(Wc),
$$

so $\Pi F=\eta F$. This remains true for nonorthogonal spatial orbitals:
orthogonalization applies the same exterior map $T$ to the two sectors.

## Hamiltonian contraction

The existing direct-CI Hamiltonian, including its one-body contribution,
is written using symmetric one-body channels as

$$
H=\mathcal E^T K\mathcal E.
$$

Here $\mathcal E$ scatters a product coefficient vector into channel values
$d_{iP}=(E_P X)_i$, with product index $i=(a,b)$ and orbital-pair channel
$P$. Each $E_P=E_{\alpha,P}+E_{\beta,P}$ commutes with exchange.
The symmetric matrix $K$ contains the two-electron integrals and absorbed
one-electron terms; it acts on channels, independently at every product point.
Thus $d_{ba,P}=\eta d_{ab,P}$.

For two vectors in the same parity sector,

$$
(Fc_L)^T H(Fc_R)
=\sum_{a,b}d_L(a,b,:)^T K d_R(a,b,:).
$$

The two members of each off-diagonal exchange orbit contribute equally:
the two signs multiply to $\eta^2=1$. Choose representatives $a\leq b$,
with weights

$$
w_{ab}=\begin{cases}1,&a=b,\\2,&a<b.\end{cases}
$$

For $\eta=-1$ the diagonal channel values vanish exactly, so only $a<b$
is needed. If $Q$ selects these representatives and $D$ stores their weights,
then the projected structure Hamiltonian is exactly

$$
H_{\mathrm{VB}}
=F^T\mathcal E^T Q^T D K Q\mathcal E F.
$$

This is a quotient of the **intermediate contraction**, not a replacement
of the structure basis by spin products. Davidson remains in structure space.
The complete one-generator support closure is first defined as before; the
quotient removes only its exchange-duplicate representatives.

## Why no mirrored output is required

The streamed scatter computes $Q\mathcal E Fc_R$. After channel multiplication
and weighting, the existing gather and structure projection evaluate

$$
F^T\mathcal E^T Q^T D K Q\mathcal E Fc_R.
$$

One must not claim that the intermediate vector before $F^T$ is the full
Hamiltonian image. It is generally not exchange symmetric. It is an exact
representative for contraction with any left vector in the prescribed sector.
The selected-structure projection interface can therefore consume this vector
directly. A full-carrier response must instead apply the sector projector
described below; the representative alone does not meet its interface contract.

Input columns needed by the original overlap action may not all be requested
by the quotient Hamiltonian sweep. Their Gram contributions must still be
evaluated. The implementation tracks consumed overlap panels separately.

## Integral adjoint and complete response

For two inputs $L$ and $R$ in the same sector, define $d_L=\mathcal E L$
and $d_R=\mathcal E R$. The adjoint of the channel kernel is

$$
\overline K=d_L^T d_R=(Qd_L)^T D(Qd_R).
$$

The existing chain rule from $K$ to the one- and two-electron integrals is
unchanged. The identity applies to accepted-state densities and to transition
densities with unequal $L$ and $R$, including orbital and CI derivatives:
spatial orbital generators act identically on both spin sectors and commute
with $\Pi$. It also applies separately to every root of a state average if
the structure basis has one common exchange parity. Mixed-parity bases and
unequal spin lists retain the unreduced contraction.

A fused response additionally requires the complete vector $HL$, not just its
projection onto the selected structures. Let

$$
P_\eta=\tfrac12(I+\eta\Pi),\qquad
z=\mathcal E^T Q^T D K Q\mathcal E L.
$$

The representative weights obey
$P_\eta Q^T DQ P_\eta=P_\eta$. Both $\mathcal E$ and $K$ preserve the
exchange sector, hence

$$
HL=P_\eta z,\qquad
(HL)_{ab}=\tfrac12(z_{ab}+\eta z_{ba}).
$$

The projection is performed in place, with one owner per unordered pair.
It needs linear carrier work and no extra full-carrier array. Omitting this
projection would give an incorrect response even though a final same-sector
scalar contraction could still pass.

The implementation packs only representative rows into each tile, including
tiles crossing the diagonal. For identical lists of size $U$, there are
$U(U+\eta)/2$ rows rather than $U^2$ per right-hand side. Prefix offsets are
arithmetic sums of row lengths; no full product index array is constructed.
This approximately halves scatter and channel-GEMM work, but does not change
their polynomial scaling. Two shared workspaces retain the same
64 MiB total bound. There is no per-thread carrier or pair-kernel adjoint.

At an accepted point $L=R$, sigma and the integral adjoint share one scatter.
The second operand is copied within the existing bounded workspace before
applying multiplicities, rather than regenerating its one-body links. Thus
the former standalone sigma plus two adjoint scatters become one link sweep.
This fusion also applies when exchange reduction is unavailable.

## Derivatives and limits

For spin-independent spatial orbital variations, $\dot T$ is identical in
the alpha and beta sectors and $\dot F$ retains parity $\eta$. Differentiating
the same bilinear identity therefore preserves the quotient for integral
adjoints and orbital derivatives; $Q,D$ are constant. For example,

$$
\delta f
=\sum_{i\in Q}w_i\left[
\delta d_{L,i}^T K d_{R,i}
+d_{L,i}^T\delta K d_{R,i}
+d_{L,i}^T K\delta d_{R,i}\right].
$$

This formula does not authorize deleting either left or right derivative.
The present implementation changes the projected H/S action only; separate
response/adjoint kernels have not yet been converted to this quotient.

The reduction is approximately a factor of two in generator target rows,
not an asymptotic removal of the spin-string product space. For a complete
product domain of dimension $U\times U$, there are $U(U+1)/2$ even-sector
representatives or $U(U-1)/2$ odd-sector representatives. The same orbit
counting applies to a symmetric restricted support closure.

Input transformation, overlap, support generation, and backprojection retain
their own costs. Their overhead can limit the wall-time benefit. No extra
full carrier or per-thread large matrix is allocated. Unequal spin sectors,
mixed parity, and arbitrary nonsymmetric expansions retain the unreduced
domain because the quotient identity does not apply to them.

## Validation contract

- Compare the parity-projected quotient image against ordinary full sigma
  for both signs, multiple RHS, and identical string lists.
- Compare selected-structure H/S images against independent occupied-minor
  reference contractions, including few structures, odd/even parity, and
  small streaming budgets.
- Retain nonsymmetric and unequal-spin tests to prevent false symmetry use.
- Check CAS(12,12) selected sets and converged 240 subspaces with unchanged
  optimizer and tolerances. Compare saved historical data, never rerun an
  old executable to manufacture a matched baseline.
