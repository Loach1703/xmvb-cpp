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

## Alpha links

For a fixed target beta group $g$, the alpha adjoint previously looked up an
original target row through the halo map, then decoded its grouped image
address separately for each excitation edge. The grouped image address is
independent of that edge and of the temporal response. Define instead

$$
J_{v g}(h)=
\begin{cases}
\pi(v N_T+i_g(a_h)), &(a_h,g)\in\mathcal T,\\
-1, &(a_h,g)\notin\mathcal T,
\end{cases}
$$

where $a_h$ is halo row $h$ and $\mathcal T$ is the current tile's actual
target set. Composing this map once gives exactly the alpha-edge contraction

$$
\sigma^{(v)}_{a g}\mathrel{+}=
\sum_{\ell\in\mathcal L_\alpha(a)}
\mathbf 1_{J_{v g}(h(a_\ell))\ge0}\,
s_\ell Z_{J_{v g}(h(a_\ell)),p_\ell}.
$$

Here $\mathcal L_\alpha(a)$ is the immutable single-replacement row of string
$a$, $a_\ell$ is its neighboring string and $h(a_\ell)$ is its halo address.
A missing halo neighbor contributes zero before querying $J$. Diagonal
occupation terms use the same composed address. Their alpha and beta
occupations and the edge signs, order and numerical additions are unchanged.
The fixed beta occupation view is obtained once per output interval rather
than inside the alpha-row loop.

The stored map now contains $wGh$ integers for $w$ RHS rather than $Gh$
original-row addresses. At $w=1$ its size is unchanged. At larger $w$ the
additional integer addresses are explicitly charged to the existing tile
budget, so the planner can reduce tile extent rather than exceed that budget.
There is still no edge-by-target table or per-worker copy. The construction
cost is $O(wN_T+wGh)$, including sentinel initialization; physical alpha-link
contributions retain their necessary work. Reusing $J$ removes the repeated
layout decode along each valid edge and across temporal responses without
changing $E^T$, the CI space, spin exchange quotient or output projection.

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

The beta schedule stores three integers per target in place of its former
one-integer target-row list. The additional two integers are charged to the
same tile admission budget, including when only one target fits. No numerical
carrier, edge-by-target table, global product mask or per-worker array is
introduced. Threads keep disjoint destination intervals and need no atomics.
