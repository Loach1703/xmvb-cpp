# Exact-demand exterior columns for tiled structure actions

## Physical operator

Let $W$ map selected VB structure coefficients to their original spin-string
products. With upper-triangular orbital factor $R$, write

$$
T_\sigma=\wedge^{m_\sigma}R,\qquad
Y=T_\alpha(WC)T_\beta^T.
$$

The selected actions are exactly

$$
HC=W^T T^T\sigma(TWC),\qquad MC=W^TT^TTWC,
\quad T=T_\beta\otimes T_\alpha.
$$

Here $W$ includes the product-index permutation for column-vectorization of
the displayed alpha-by-beta matrices. $H$ and $M$ act only in the selected
structure space; $\sigma$ is the
orthonormal direct-CI Hamiltonian action. The formulas do not assume that
the selected structures are closed under the exterior map or Hamiltonian.
The current change leaves these operators and their derivatives unchanged.

## Requested columns, not enclosing panels

For the exact set $B$ of beta columns requested by one sigma halo, and the
original touched products $g=(a_g,b_g)$, expansion can be evaluated as

$$
V_{a,j}=\sum_{g:a_g=a}(WC)_g(T_\beta)_{j,b_g},\qquad
Y_{:,B}=T_\alpha V_{:,B},\qquad j\in B.
$$

Neither expression uses beta rows outside $B$. Generating their enclosing
contiguous panels is unnecessary. A request may be split into bounded batches
of missing columns without changing any sum. Hits are read before insertion,
so later evictions cannot invalidate borrowed input pointers. Cache residency
is charged per column; generation and parallelism remain batched.

For an output panel $J$, its contribution to the touched-product image is

$$
h_g\mathrel{+}=\sum_{j\in J}
 (T_\alpha^T Z_{:,J})_{a_g,j}(T_\beta)_{j,b_g}.
$$

This is linear in the partial sigma output $Z$. Eviction may therefore project
partial outputs exactly; no physical term is omitted. Output residency is
also per actually written beta column. Evicted columns are moved, not copied,
into a bounded pending batch and projected together. Unwritten columns of
their former enclosing panels are never allocated, transformed or projected.
Real partial columns may still require more than one projection if revisited;
this change does not claim to eliminate that time/storage tradeoff.

The sigma writer borrows a budget-fitting group of most-recently-touched
column owners before launching one shared task group. Their total charged
payload fits the LRU, so inserting later columns in that group cannot evict
an earlier borrowed owner. Workers then update disjoint beta/RHS/alpha ranges.
This removes per-column OpenMP joins without copying the coefficient arrays.

The overlap image uses $Z=Y$. Each input column contributes to overlap exactly
once per RHS block, even if it is later regenerated. A mixed batch masks the
already-consumed columns before its overlap projection. This is an exact
partition of the contraction, not an amplitude cutoff or stationarity test.

## Direct exterior-row recursion

For occupations $I=(i_1<\cdots<i_m)$ and $J=(j_1<\cdots<j_m)$,

$$
T_{IJ}=\det R[I,J].
$$

The former row builder propagated a unit row through every elementary shear
of the complete single-spin space. Instead, prepend occupied rows in
descending order. Set $D_0(\varnothing)=1$ and define

$$
D_k(J)=\det R[(i_{m-k+1},\ldots,i_m),J].
$$

Laplace expansion along the newly prepended row gives

$$
D_k(j_1,\ldots,j_k)=
\sum_{r=1}^{k}(-1)^{r-1}R_{i_{m-k+1},j_r}
D_{k-1}(J\setminus\{j_r\}).
$$

Upper triangularity implies a strict structural zero unless

$$
j_r\ge i_{m-k+r},\qquad r=1,\ldots,k.
$$

Only combinations satisfying this condition are enumerated. Since
$i_{m-k+1}\ge m-k$ in zero-based indexing, stage $k$ only needs the
$k$-subsets of sites $m-k,\ldots,n-1$. The stage capacity is

$$
\binom{n-m+k}{k}\le\binom{n}{m}=U_\sigma.
$$

Two single-spin arrays are sufficient even above half filling. No inverse,
deleted-minor factorization, dense compound matrix, or product-sized table is
introduced. Empty and fully occupied spin sectors obey the same recurrence.

The implementation ranks combinations in colexicographic order,

$$
\operatorname{rank}(j_0,\ldots,j_{k-1};b)
=\sum_{r=0}^{k-1}\binom{j_r-b}{r+1},\qquad b=m-k.
$$

Predecessor ranks after removing each occupied site are obtained from prefix
and suffix sums of Pascal coefficients, taking $O(k)$ work for all removals,
not $O(k^2)$. The caller's actual spin-string ordering is restored through
a linear-size rank permutation. Both stage arrays are cleared only over the
current stage capacity so invalid predecessor states are exact zeros.

For one requested row $I$, if $N_k(I)$ combinations pass the structural
condition, arithmetic work is

$$
O\!\left(\sum_{k=1}^{m}kN_k(I)\right),
$$

versus scanning up to $\binom{n}{2}\binom{n-2}{m-1}$ shear edges per
row in the previous builder. A worst-case upper bound on the recurrence
work is

$$
\sum_{k=1}^{m}k\binom{n-m+k}{k}
=(n-m+1)\binom{n+1}{m-1}.
$$

This is a row-generation reduction, not a claim to remove the many-electron
exponential dimension. Actual input regeneration also depends on the sigma
access graph. Removing unrequested columns and reducing miss cost are separate
improvements.

## Supported shear graph and RHS batching

Each shear edge $(s,t)$ is relevant to a beta column only if both ends belong
to that column's exact allowed alpha set $A$. Store the shear's source set
$P$ as a shared bitset. Visit set bits of $A\cap P$; the edge index is its
source-bit rank in $P$, using one stored prefix count per word and a
population count within that word. The target is then tested against $A$.
This is precisely the original predicate without scanning every rejected
edge. Metadata are $O(n^2U_\alpha/64)$ words, shared by all columns and threads.

Different RHSs have the same allowed graph. Apply each relevant edge to all
RHS lanes before advancing to the next edge. Shear/edge ordering within each
individual RHS remains unchanged; only repeated graph traversal is removed.

Full-product residency is not required for batching RHSs. The existing
channel-workspace limit determines RHS width; the full-product residency
estimate is used only by the already-resident contraction order. The sigma
executor still derives its actual halo/channel batch size from the unchanged
workspace bound.

## Storage and correctness contract

The Pascal table is $O(n^2)$ and the ordering map is $O(U_\sigma)$.
Worker scratch is two $U_\sigma$ arrays in one operation-owned allocation;
the worker count is bounded by
the existing 64 MiB scratch budget, as well as available threads and requested
rows. Workers never copy integrals, CI vectors, or a spin-product matrix.
If even two single-spin vectors exceed the scratch target, one worker is
the unavoidable minimum; no fixed budget is used to reject the physical
representation. For C16 at 32 workers, these slices total about 6.28 MiB.
Generation batches use the existing panel extent and workspace budget;
column-resident primal data still have the same 64 MiB cap.

No cache budget, convergence threshold, selected space, preconditioner or
CIAH strategy is changed. Correctness is checked against independent compound
matrix minors, ordered-shear actions, complete H/S structure actions,
arbitrary selections, unequal spins, and existing derivative/transport tests.
The tiny-residency reader test forces evictions and checks that only requested,
nonzero-support columns are generated, including multiple RHSs and execution
inside an existing OpenMP team.

Accepted-point jet residency is not silently declared solved. Its separate
cost must be measured after this initial-CI repair; the jet still has its own
panel sweep. This work removes enclosing-panel replay in structure H/S actions,
not every possible regeneration in all gradient/HVP consumers.
