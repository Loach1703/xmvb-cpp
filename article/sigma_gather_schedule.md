# Exact composition of the tiled sigma adjoint schedule

The orthogonal spin-string Hamiltonian is applied as $E^T K E$. The selected
VB coefficients and Davidson vectors remain in the original structure basis;
this section changes only the evaluation schedule of the generator adjoint
$E^T$. It assumes neither closure nor completeness of the selected structures.

## Beta links

For one beta excitation edge $\ell:g\to b$, its sign $s_\ell$ and packed
orbital-pair channel $p_\ell$ do not depend on alpha. Let $T_g$ be the actual
alpha targets of group $g$ in the current tile, and let $D_b$ be the caller's
output support. Its contribution to RHS $v$ is exactly

$$
\sigma^{(v)}_{a b}\mathrel{+}=
s_\ell Z_{\pi(v N_T+i_g(a)),p_\ell},
\qquad a\in T_g\cap D_b.
$$

Here $N_T$ is the total target count in this tile, $i_g(a)$ is the original
tile-row address of target $(a,g)$, $Z$ contains the channel images, and $\pi$
is their channel-group layout permutation. Different beta edges generally
have different $T_g\cap D_b$; dropping those intersections is not an identity.

The layout groups all image rows by source-channel envelope. If target $i$
has envelope $c(i)$ and that envelope contains $m_{c(i)}$ targets per RHS,
the existing layout has the exact affine relation

$$
\pi(v N_T+i)=\pi(i)+v m_{c(i)}.
$$

Consequently each target needs only one symbolic record

$$
\bigl(a,\pi(i_g(a)),m_{c(i_g(a))}\bigr).
$$

These records are sorted by alpha within each beta group. Sorting already
ordered columns is unnecessary. Each incoming beta edge intersects this
sorted list directly with the destination rows and evaluates the affine
address, without a halo-index lookup, a target-index lookup and another
layout-permutation lookup per contribution. The RHS stride is the envelope
population, **not** the beta-group population. With one RHS it can be zero.
Every temporal response uses the same records but retains its own channel
values and cotangents.

## Strictly removed work and storage bound

The generator already supplies every target. Constructing a sorted target
list by scanning all $h$ halo rows for each of $G$ beta groups is therefore
unnecessary. The schedule now visits the $N_T$ supplied targets, sorting only
an unordered group. Already sorted columns take $O(N_T)$ preparation rather
than $O(Gh)$; arbitrary unordered columns additionally require their actual
group sorts, not an assumption of covalent structure.

The address composition is performed once per tile, instead of once per
participating beta edge, RHS and temporal response. This does **not** remove
the physically necessary edge contributions or change their accumulation
order, channel multiplication, overlap action, orbital response or optimizer.
It is a reduction in repeated sparse address work, not a claimed lower
universal power of the active-space size.

The current tile stores three integers per target in place of its former
one-integer target-row list. The additional two integers are charged to the
same tile admission budget, including when only one target fits. No numerical
carrier, edge-by-target table, global product mask or per-worker array is
introduced. Threads keep disjoint destination intervals and need no atomics.
