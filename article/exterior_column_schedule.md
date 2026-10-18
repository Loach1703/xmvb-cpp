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

## Last-consumer scheduling

The intermediate channel target $(a,b)$ reads and writes beta columns in
$N[b]$, the identity plus the one-replacement neighbors of $b$. For the
actual beta-target order $b_0,\ldots,b_{L-1}$, define

$$
\ell(j)=\max\{t:j\in N[b_t]\}.
$$

This is an integer lifetime bound derived from the generator graph, not an
amplitude threshold. Computing it requires one pass over the single-spin
links and $O(U_\beta)$ integers; neither channel amplitudes nor a product
table are stored. It remains a valid conservative bound for arbitrary
selected support, unequal spins and exchange-representative sweeps.

After the final alpha chunk of target $b_t$ has been consumed, a column
with $\ell(j)\le t$ cannot participate in any remaining scatter or adjoint
scatter. Release its input immediately. Consume and release its output;
its last partial contribution can join the existing bounded projection
batch. If several alpha chunks belong to one target, the target remains
live until the cursor is exhausted. Parallel workers and writer tasks must
finish before release.

This eliminates the mathematically unnecessary residency of dead columns.
Previously a final-use column could become most recently used in the LRU,
displacing a live column whose transform then had to be repeated. Removing
a dead entry cannot introduce a later input miss. Output release does not
add a projection: the same completed partial image otherwise had to be
projected on eviction or at the final flush. This does **not** eliminate
pressure-driven eviction of still-live columns; the latter requires a
different scheduling argument and must be measured separately.

## Share orbital factors, not state-dependent columns

The same exterior row $t_j=(T_\beta)_{j,:}$ occurs in both input construction
and the projected output sum above. It is independent of $C$, the channel
kernel, and the Davidson RHS. Its value therefore remains identical across
both consumers and every CI iteration at a fixed orbital point.

The pre-existing 64 MiB exterior-row cache now admits exact requested
single rows under one shared key. Missing rows are still generated as a
parallel batch; cache admission remains single-row. Contiguous row-block
consumers share the same total budget, not a second factor cache. Returned
matrices are operation-local copies, so eviction cannot invalidate sigma
views. A new orbital transform creates a new cache, preventing reuse across
different orbital points. This removes duplicate factor generation on hits
without storing a complete exterior matrix or a spin-product grid. The
state-dependent input and output column budgets are unchanged.

Generated minors, cached single rows, and intermediate row assembly now all
use row-contiguous storage. Only the completed requested-row panel is
converted to the column-major view consumed by the source contraction. This
removes an intermediate full-panel conversion and repeated strided row
stores. It changes neither a minor recurrence nor a contraction sum.

## Contract only nonzero beta factors

For a requested beta-row set $B$, define

$$
J_b=\{j\in B:(T_\beta)_{j,b}\ne0\}.
$$

The input scatter and output contraction use only these entries:

$$
V_{a_g,j}\mathrel{+}=(WC)_g(T_\beta)_{j,b_g},\qquad j\in J_{b_g},
$$

$$
h_g\mathrel{+}=\sum_{j\in J_{b_g}}
(T_\alpha^T Z_{:,B})_{a_g,j}(T_\beta)_{j,b_g}.
$$

Each bounded sampled-factor panel is stored in compressed sparse column
(CSC) format, indexed by the original beta string $b$. Two parallel passes
count and pack disjoint CSC columns; only their offset prefix sum is serial.
The conversion removes only values equal to floating-point zero, not values
below a magnitude tolerance. Tiny nonzero values, including subnormals, are kept.
The same CSC panel serves every RHS and the input overlap consumer; there is
no allocation per touched product.

For an upper-triangular $R$, $T_{I,J}=0$ unless $j_r\ge i_r$ for every
occupation. In the complete $m$-electron spin sector the number of symbolically
allowed pairs is

$$
N_{\mathrm{nz}}^{\mathrm{sym}}(T)
=\binom{n}{m}^2-\binom{n}{m-1}\binom{n}{m+1}
=U_\sigma^2\frac{n+1}{(m+1)(n-m+1)}.
$$

The subtraction follows by reflecting the two occupation-prefix paths at
their first ordering violation. At $(n,m)=(16,8)$ the full-spin factor density
is at most $17/81$, about 21.0%; actual sampled panels can have a different
density and additional exact zeros. This counting does not assume closure of
the selected VB structures.

If $G$ is the touched original-product set, the two contractions formerly
visit $O(b|G||B|)$ entries. CSC consumption visits
$O(b\sum_{g\in G}|J_{b_g}|)$ entries instead. The current conversion still
scans a bounded dense factor panel in $O(U_\beta|B|)$, and the row recurrence
still has its stated generation cost. Thus this is not a claim that the
entire H/S action has acquired the sparse-consumption complexity. Storage
for CSC is proportional to its panel nonzeros plus $U_\beta$ offsets; neither
a complete exterior matrix nor a complete product grid is retained.

## One-pass target planning and a shared resource budget

Let $I_t$ be the exact generated target set in a tile, $A_t$ its distinct
alpha strings, and $B_t$ its distinct beta strings. With $N_\sigma(i)$ denoting
the identity plus single-replacement neighbors, the input rectangle is

$$
\mathcal H_t=\bigcup_{a\in A_t}N_\alpha(a),\qquad
\mathcal D_t=\bigcup_{b\in B_t}N_\beta(b).
$$

These are unions, not the product of a target count and a worst-case degree.
A persistent support cursor emits each target once. A trial target adds
only previously unseen neighbors; if admitting it would exceed the numeric
payload budget, that one target is held for the next tile. The cursor is
never rewound. Halo sorting, coordinate-map construction and target packing
are performed once per completed tile, rather than in separate repeated
planning passes. This changes the partition, not the target set, its
support labels, or exchange multiplicities.

For $r=|I_t|$, $h=|\mathcal H_t|$, $d=|\mathcal D_t|$, $g=|B_t|$,
RHS width $b$, channel count $P$, and $s$ simultaneously live channel slabs,
the planned numeric payload is

$$
B_t^{\mathrm{num}}=8\left(bhd+sbrP\right)+4hg.
$$

The final term is the target-index map. For an in-place temporal response,
the input allocation is instead $8\max(bhd,\min(n_{\mathrm{workers}},br)P)$;
the already-consumed halo supplies the GEMM scratch. No channel slab is
replicated by worker count. Single-spin links, cursors, and small metadata
are separate from this numeric bound.

The H/S action now defaults to a **4 GiB process-shared tile budget**,
configurable with `--ci-workspace-mib`. Input and output column residency
receive one quarter each; sigma receives one quarter; the remaining quarter
is reserved for transformation panels. Allocation remains demand-driven:
small carriers do not allocate the full allowance. Exterior-factor caches
are immutable accepted-point data with their existing separate bounded
residency; AO integrals, Davidson vectors and library/allocator storage are
also outside this tile-payload budget. It is therefore not a process RSS
limit. The resource option changes neither precision nor the variational
space. The product-jet derivative consumers use the same configured allowance
with the partition below, not an independent allowance for each consumer.

Column retirement still uses exact last generator use. Its sorted order is
built once and advanced monotonically, avoiding a scan of all beta strings
at every tile boundary. The numerical workspace retains one shared idle
allocation, never one allocation per OpenMP worker; a smaller-budget caller
discards an oversized idle allocation before beginning its action.

## Completed images before reverse contraction

Write a normalized selected state as $y=TWc$. Its orbital variation obeys

$$
\mathrm d E=y^T(\mathrm dH)y+2[(H-E)y]^T(\mathrm dT)Wc.
$$

CI stationarity gives $W^TT^T(H-E)y=0$, not $(H-E)y=0$ in an arbitrary
selected space. The second term therefore remains in both gradient and HVP.
The scheduling change below uses linearity, not a small-residual approximation.

For an output column panel $b$, let $z_{bt}$ be the image contribution from
generated sigma tile $t$, and let $\mathcal G(z,y)$ be the generator cotangent
with its right argument held fixed. Then

$$
z_b=\sum_t z_{bt},\qquad
\sum_t\mathcal G(z_{bt},y)=\mathcal G(z_b,y).
$$

The same identity applies to exterior-adjoint projection and to each term of
its directional derivative. Repeated consumption of partial output panels
is mathematically correct, but needlessly repeats neighbor reads, exterior
generation and projection. Prefer a completed image and consume it once.

Sigma retirement now extends to adjoint, directional and temporal-response
sweeps: a column retires only after every response and forcing writer in its
last tile has finished. A panel retires after all its columns. Once a sweep
has capacity evictions, completed output panels are immediately consumed and
released, reducing pressure on unfinished panels. Fully resident sweeps keep
their completed panels until exchange projection and accepted-state retention.
An input panel must NOT be retired using this shorter lifetime: generator
adjoints still read neighboring inputs while consuming completed outputs.
Retirement runs inside the final writer's existing OpenMP single region.
Testing a mutable retirement list outside it can race the next tile planner
and send workers to different conditional barriers; no additional barrier
is needed when retirement shares the final writer's synchronization.

With total configured budget $B$, at most $B/4$ is reserved for accepted
panels across all roots. Divide the remainder into four equal parts for
input panels, output panels, sigma, and reverse-transform workspace. In the
last part, projection accumulation and directional exterior rows each have
one quarter; bounded transient panels have the other half. Small explicit
panel caps can reduce residency but never alter support or derivatives.
The panel extent respects the transient bound as well as spin-space sizes.

Accepted $y$ and restricted $(H-E)y$ share one structural index pattern.
For $N$ supported coordinates, their payload is approximately $20N$ bytes,
not twice the one-lane $12N$ bytes. Retention must also fit the input-panel
allowance for a subsequent two-lane coupled direction; otherwise borrowing
all directional panels across a parallel reader would invalidate LRU views.
For the C16 support $N=34763300$, the payload excluding offsets is663.06 MiB.
No complete channel array is retained. Immutable exterior row caches,
integrals, spin links, metadata and allocator residency remain outside this
numeric tile accounting and must be included in measured process RSS.

For states larger than the allowed live frontier, exact partial-image
consumption is still required. This change does not claim that an arbitrarily
small budget can avoid all recomputation, or that a sigma tile budget is a
bound on whole-process RSS.

## Restricting channels without dropping exterior images

For a generated product node $t$, define the channel vector
$a_t=(E y)_t$, where $E$ collects the packed spin-free one-body generators.
The two-body channel operation and its integral adjoint are

$$
c_t=K a_t,\qquad G=\sum_t w_t a_t b_t^T.
$$

Here $b_t=(E z)_t$ is the channel vector of the other argument and $w_t$ is
the exchange-orbit weight. If a Boolean selector $J_t$ contains every
structurally reachable input channel, then

$$
a_t=J_t\widetilde a_t,\qquad
c_t=(KJ_t)\widetilde a_t,\qquad
G=\sum_t w_tJ_t\widetilde a_t\widetilde b_t^TJ_t^T.
$$

The second identity retains **all output channels**, including those leading
outside the selected input support. It does not assume closure and does not
set the exterior residual to zero. An output selector needs its own generator
connectivity proof; it cannot simply be copied from the input selector.

For a complete covalent seed family, let

$$
d_k(t)=N_\alpha(t;0{:}k)+N_\beta(t;0{:}k)-(k+1).
$$

Depth-zero support has every $d_k\geq0$. At a depth-one frontier node outside
that support, the minimum is $-1$. Let $f$ and $l$ be the first and last
negative prefixes. Every single replacement returning to depth-zero support
must raise an electron from $q$ to $p$ with

$$
0\leq p\leq f\leq l<q<n,\qquad k_t=(f+1)(n-1-l).
$$

This rectangle is a safe channel envelope; occupations can remove additional
channels. Supported nodes retain all $P=n(n+1)/2$ channels. Grouping nodes
with the same $(f,l)$ allows an ordinary dense GEMM with $KJ_t$, changing its
channel arithmetic from $O(m_tP^2)$ to $O(m_tPk_t)$. This is an exact
support-dependent work reduction, not a smaller asymptotic power in $n$ for
every family. The potentially smaller $O(m_tk_t^2)$ adjoint is a subsequent
step, not part of the initial implementation.

The current envelope certificate is specific to a complete covalent seed
family. Arbitrarily selected and ionic seed families keep their original
exact channel operation until their own channel envelopes are constructed.
The grouped multiplication interface itself is independent of the family.
Kernel panels are shared, charged against the existing action budget and
never replicated per worker. Packing reuses disjoint pieces of the consumed
input halo after a team barrier; source and output channel slabs are unchanged.

An exhaustive C8 diagnostic and a deterministic C16 stride-257 census compare
the envelope with explicitly enumerated supported one-generator neighbors.
For C16,197498 sampled nodes out of50756849 exchange-quotient targets have
mean exact input-channel count39.2394, versus136 unrestricted channels.
The conservative envelope mean is60.7388. Its predicted channel multiplication
work is44.66% of unrestricted work. This is a FLOP estimate, not a measured
sigma or SCF speedup; packing, support traversal and exterior transforms remain.

## Ragged channel storage and fused directional images

The first grouped-GEMM implementation still built dense channel rows and
packed/unpacked them around each multiplication. Its C16 one-step runtime
did not improve. The next implementation places the rows of each selector
group directly in contiguous, variable-width matrices:

$$
A_g,B_g\in\mathbb R^{m_g\times k_g},\qquad
C_g=A_g(KJ_g)^T\in\mathbb R^{m_g\times P}.
$$

The scatter addresses $A_g$ directly, skipping only generator edges whose
channels are certified absent by $J_g$. There is no dense intermediate input,
no numerical row permutation and no post-GEMM unpacking. Integer row maps
connect the original product enumeration to the grouped output. All output
channels remain present, so the existing complete-image adjoint is retained.

The reverse contraction now also uses the smaller matrices:

$$
G=\sum_g J_g\left(A_g^TW_gB_g\right)J_g^T,
$$

where $W_g$ contains exchange-orbit multiplicities (one on the diagonal,
two otherwise). A shared accepted-state channel matrix remains unweighted;
weights act only on the Gram and the completed sigma image. This permits
its reuse across temporal responses and directional Hamiltonians. When the
bra and ket coincide, the Gram borrows the same channel storage twice rather
than copying or scattering it again.

For a directional image, perform the second product as an accumulation:

$$
\delta C_g=A_g(KJ_g)^T+B_g(\delta KJ_g)^T.
$$

Unless a consumer explicitly requests the independent forcing image,
the second term needs no full-size temporary. An explicit forcing request
reuses the consumed source storage only after every worker has completed
the base product. Grouped BLAS and original-row scatter have different row
owners; their synchronization is explicit, not inferred from equal row counts.

With $R$ target/RHS rows and $L=\sum_gm_gk_g$, the simultaneous channel
payload is $L+RP$ doubles for H/S or a single adjoint, $2L+RP$ for a shared-state
direction/response sequence, and $L+2RP$ when a separate forcing image is
requested. The planner computes this payload while generating nodes, also
charging the row maps and existing target map. No per-thread numerical copy
of these slabs is introduced. Small per-worker Grams borrow disjoint portions
of the consumed halo; if the halo cannot hold them, workers own independent
global Gram entries without allocating extra scratch.

The unrestricted family is the same algorithm with $J_g=I$ and implicit
identity row maps. The earlier full-channel packing implementation and its
in-place expansion logic are removed, rather than kept as a second solver.
