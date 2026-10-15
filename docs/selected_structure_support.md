# Exact orthogonal-CI action for selected VB structures

## Variational space and computational space

Let $P$ map coefficients of the selected VB structures to their original
alpha/beta spin-string products. Its nonzero rows form a set $G$. Adding
indices for other spin strings does not add columns to $P$ and therefore
does not enlarge the variational space. With $\Phi=\chi R$ and upper-triangular
$R$, define $T_\sigma=\wedge^{n_\sigma}R$. For a structure vector $c$,

$$
X=\operatorname{reshape}(Pc),\qquad
Y=T_\alpha X T_\beta^T.
$$

Writing $F:c\mapsto Y$, the required actions are

$$
H_{\mathrm{VB}}c=F^\dagger H_{\mathrm{orth}}Fc,
\qquad Mc=F^\dagger Fc.
$$

No invariance or closure of the selected structure space under $T_\sigma$ is
assumed. In particular, intermediate spin products outside $G$ are required.
The Davidson vectors remain in the selected structure space throughout.

## Exact structural support

For a spin string $A$, let $p_k(A)$ count its occupied orbitals among indices
$0,\ldots,k-1$. Define the downward ideal

$$
D(A)=\{a:p_k(a)\ge p_k(A)\ \text{for every }k\}.
$$

Equivalently, the ordered occupied indices satisfy $a_i\le A_i$. Every
nonzero upper-triangular minor respects this order. Thus an exact structural
support for the transformed selected space is

$$
\mathcal S=\bigcup_{(A,B)\in G}D(A)\times D(B).
$$

Here "structural" permits coefficients to cancel, and permits particular
orbital transformations to contain additional zeros. Such numerical zeros
are not used for screening. Membership is exact for the generic triangular
transform of the original product support. It is sufficient for any selected
linear combinations on that support, including covalent and ionic structures.

Both spin conditions must refer to the **same** seed $(A,B)$. Taking the
Cartesian product of the separate marginal unions generally overestimates
the reachable set. It is safe but loses the intended computational reduction.

For the complete singly occupied seed family, the union simplifies to

$$
\mathcal S=\{(a,b):p_k(a)+p_k(b)\ge k\ \text{for every }k\}.
$$

That closed form need not be the exact support of a selected subset.

## Intermediate space of the Hamiltonian

With $E_{pq}=a^\dagger_{p\alpha}a_{q\alpha}
+a^\dagger_{p\beta}a_{q\beta}$, the spin-free Hamiltonian contains products
of two one-body generators. The first generator can leave $\mathcal S$ even
when the second returns to it. Define

$$
D_1(A)=\{a:p_k(a)+1\ge p_k(A)\ \text{for every }k\}.
$$

This is exactly the union of $D(A)$ and its single-excitation neighbors. One
replacement can change any prefix count by at most one. Conversely, if all
prefix deficits are at most one, move the occupied site immediately after
the last negative-prefix interval into the empty site at the beginning of
the first such interval. The resulting string belongs to $D(A)$.

Consequently the required intermediate support is

$$
\mathcal I=\bigcup_{(A,B)\in G}
\left[D_1(A)\times D(B)\;\cup\;D(A)\times D_1(B)\right].
$$

Only one spin sector may use the depth-one allowance, since a single
generator acts on one spin sector. Allowing one replacement in both sectors
would describe a larger, unnecessary intermediate space.

For each generator channel, form $Z_{pq}=E_{pq}Y$ only on $\mathcal I$.
Contract the integral kernel between channels there, then apply the adjoint
generators and retain only $\mathcal S$. Thus the projected Hamiltonian
action retains all two-generator paths and is exactly $P_{\mathcal S}
H_{\mathrm{orth}}P_{\mathcal S}$, including paths that leave $\mathcal S$.

## Streaming transforms and projection

Every elementary forward triangular shear moves an occupied orbital to a
lower index. The set $\mathcal S$ is closed under these moves, so all forward
paths from $G$ stay within $\mathcal S$. The adjoint reverses those paths.
Every reverse path contributing to $G$ is the reversal of a forward path;
discarding other adjoint paths does not change $F^\dagger$ on the selected
structures. Tiles may therefore restrict both transforms to $\mathcal S$.

The implementation should expand $Pc$ once per RHS batch, consume tiles in
$\mathcal S$ and $\mathcal I$, accumulate on the touched original products,
and apply $P^T$ once. It must not form a dense matrix of structural masks or
re-expand every selected structure independently for each tile.

## Direct generation, rather than product-space membership tests

Let $m$ be the number of active spatial orbitals and $K=|G|$ the number of
distinct touched original spin products. $K$ is not the number of VB structures.
All contributions from selected structures to the same original spin product
are first combined through $Pc$. Support generation must not repeat this
expansion for each structure or for each integration tile.

Write a product as a sequence of local occupation pairs,

$$
z=((a_1,b_1),\ldots,(a_m,b_m)),\qquad
(a_k,b_k)\in\{0,1\}^2,
$$

with fixed sums $n_\alpha$ and $n_\beta$. After $k$ sites, let
$u_k=\sum_{j\le k}a_j$ and $v_k=\sum_{j\le k}b_j$. A partial path carries
sets of still-compatible **paired** seeds:

$$
Q_k^{(d_\alpha,d_\beta)}=
\left\{(A,B)\in G:
u_j+d_\alpha\ge p_j(A),\quad
v_j+d_\beta\ge p_j(B),\quad 1\le j\le k\right\}.
$$

Use $(d_\alpha,d_\beta)=(0,0)$ for $\mathcal S$, and the union of the modes
$(1,0)$ and $(0,1)$ for $\mathcal I$. Keep $(0,0)$ as an additional label
when traversing $\mathcal I$. The allowance is fixed along a path: changing
which spin receives the allowance at successive sites would be incorrect.

Precompute the seed sets

$$
L_{\sigma,k}(t)=\{g\in G:p_k(g_\sigma)\le t\}.
$$

Each of the four possible local extensions updates the current set by

$$
Q_{k+1}^{(d_\alpha,d_\beta)}=
Q_k^{(d_\alpha,d_\beta)}
\cap L_{\alpha,k+1}(u_k+a_{k+1}+d_\alpha)
\cap L_{\beta,k+1}(v_k+b_{k+1}+d_\beta).
$$

A branch is absent if all its compatible-seed sets are empty, or if its
remaining sites cannot supply the required spin electron counts:

$$
0\le n_\alpha-u_k\le m-k,\qquad
0\le n_\beta-v_k\le m-k.
$$

These are exact integer conditions, not coefficient screening or a choice
between physical algorithms. The result is a prefix tree of **products**,
not separate trees independently expanded for every seed. A leaf has a unique
occupation sequence and is emitted exactly once, even if thousands of seeds
reach it. There is no post-generation hash table of all products and no
$U_\alpha U_\beta$ mask.

### Why pruning is complete

Every reachable product satisfies the constraints at each prefix, so its
path cannot be pruned. Conversely, a surviving complete path has a compatible
seed and hence belongs to the required union of downward ideals. Moreover,
each surviving partial path can be completed: place each spin's remaining
electrons as early as possible. For a surviving seed, this maximizes every
future prefix count and satisfies its remaining lower bounds. This argument
uses the complete fixed-electron-number spin universe, not a carrier with
arbitrarily deleted strings.

Thus the tree contains no surviving branches with no output below them.
If $N_d$ products are emitted at depth allowance $d$, the number of surviving
prefixes is at most $(m+1)N_d$. At most four children are tested per prefix.
This is output-sensitive enumeration, rather than a scan of the complete
Cartesian product.

### Representation and exact complexity

The compatible-seed sets are bitsets over $K$ seeds, not over spin products.
Let $W=\lceil K/w\rceil$, with machine-word width $w=64$ in the prototype.
The tables $L$ can be constructed by prefix accumulation over the count $t$.

| Operation | Time | Additional storage |
| --- | --- | --- |
| Deduplicate original seed pairs | $O(K\log K)$ | $O(K)$ |
| Build threshold bitsets | $O(mK+m^2W)$ | $O(m^2W)$ words |
| Generate $\mathcal S$ or $\mathcal I$ | $O(mWN_d)$ worst case | $O(mW)$ words for the traversal |
| Emit a bounded tile | $O(B)$ callback/storage operations | $O(B)$ product records |

The factor $W$ is explicit: this construction removes full-product scanning,
but does not make seed bookkeeping free. The actual work is proportional to
the number of visited prefixes and live-bitset word operations. No complete
single-spin index tables are required to generate occupation masks; a later
consumer may use its existing single-spin index maps.

The full covalent-family prefix rule is an algebraic simplification of this
same construction and removes the seed-bitset factor altogether. It is not
a memory-dependent change of physical representation. Arbitrary selected
families retain the general paired-seed construction.

## From generated support to numerical contraction paths

Support generation alone is not a Hamiltonian algorithm. In particular,
generating products independently for each seed and subsequently applying
the Hamiltonian would repeat physical work and is not the intended method.

### Orbital transformation

Factor each triangular exterior transform into elementary shears and a
diagonal scaling. An elementary shear contributes at most one signed
single-spin excitation from each source product. In the time-ordered graph
of shear stages, accumulate all incoming coefficients at a node before
propagating them to the next stage. Do not enumerate every complete shear
path separately. Since $\mathcal S$ is downward closed, every contributing
forward path lies inside it. Reversing this identical graph gives its exact
adjoint and hence the projection back to the original selected structures.

### Hamiltonian

For the ordered spin-free generators $E_{pq}$, define

$$
G_{pq}(i,s)=\langle i|E_{pq}|s\rangle,
\quad i\in\mathcal I,\quad s\in\mathcal S.
$$

These links follow directly from occupied-to-unoccupied single-spin
substitutions, plus diagonal occupation links. Each generator changes only
one spin string; no Cartesian product of two link lists is required. The
Hamiltonian identity is

$$
\hat H_{\mathrm{act}}=
\sum_{pq}\bar h_{pq}E_{pq}
+\frac12\sum_{pqrs}(pq|rs)E_{pq}E_{rs},\qquad
\bar h_{ps}=h_{ps}-\frac12\sum_q(pq|qs).
$$

For real left and right carriers $L,R$ supported on $\mathcal S$, set
$Z^R_{pq}=G_{pq}R$ and $Z^L_{pq}=G_{pq}L$. Then

$$
L^T H_{\mathrm{act}}R=
\sum_{s\in\mathcal S}\sum_{pq}\bar h_{pq}L_s Z^R_{pq}(s)
+\frac12\sum_{i\in\mathcal I}\sum_{pqrs}
(pq|rs)Z^L_{qp}(i)Z^R_{rs}(i).
$$

The transpose $qp$ on the left is required because $E_{pq}^\dagger=E_{qp}$.
This formula can be evaluated one generated intermediate tile at a time.
The equivalent sigma uses the same generator links in reverse. The code's
symmetric packed channels absorb the one-body term into the pair kernel;
the ordered formula above is useful as an independent reference.

**Deduplicate products, not physical contributions.** Distinct sources or
generator channels reaching one intermediate must still be added. Only the
integral contraction on the accumulated intermediate is performed once.
Dropping intermediates outside $\mathcal S$ is not valid.

### Contraction scaling and residency

Let $s=|\mathcal S|$, $i=|\mathcal I|$, $b$ be the RHS count, $d=O(m^2)$ the
number of generator channels, $e_T$ the retained shear edges, and $e_G$ the
retained generator edges. Excluding regeneration, the work has the form

$$
O\!\left(
b\,\mathrm{nnz}(P)+b e_T+b e_G+b i d^2
\right),
\qquad e_T=O(m^2s),\quad e_G=O(m^2s).
$$

The dense exact-integral channel contraction remains $O(b i m^4)$ in this
representation. Support generation reduces the carrier factor from the full
product to the generated intermediate set; it does not reduce the integral
power of $m$. A factorized integral kernel changes the channel contraction
independently, not the support theorem.

A numerical intermediate tile needs $O(Bdb)$ scalars, rather than $O(idb)$.
This alone does not bound the full pipeline: input transform panels, reverse
accumulations, link metadata, and their lifetimes must be scheduled together.
If limited residency regenerates a transform panel or repeats support
traversals, the corresponding replay factors multiply those work terms.
They must be counted, not hidden behind the word "tile". Storing the full
prefix tree or all excitation edges is unnecessary; traversal and edge
generation can be streamed.

## Limits and qualification

These formulas hold for any number and choice of selected VB structures;
Davidson remains in their coefficient space. They do not require selected
structures to be closed under orbital transformations.

The support is exact for independently variable coefficients on $G$. Fixed
relations imposed by $P$, spin symmetry, or special orbital matrices can
create additional exact cancellations. The union is then a sufficient
support, not a proof that every retained entry is physically nonzero.

There is no universal sub-product-size guarantee. A single seed with both
spin strings occupying their highest available orbital indices has a full
downward ideal in each spin sector. This is already a counterexample to a
complexity bound based only on the number of structures. Additional tensor
compression or analytically established cancellations would be a separate
algorithm, not a consequence of triangular support alone.

The streaming prefix generator is implemented in
`src/mcscf/ci/product_paths.hpp`. Small independent enumeration checks exact
support, one-generator closure, unequal spins, ionic products, shared seeds,
and empty/complete cases. An ordered-generator Hamiltonian contraction checks
that paths outside the output support are retained correctly. CAS(12,12)
tests record generated counts, time, metadata and traversal memory. The
ordered and packed contractions share the one-body-generator implementation;
this checks support and contraction consistency, not an independent fermionic
sign implementation. The test also requires a nonzero Hamiltonian contribution
from $\mathcal I\setminus\mathcal S$, so omission of the halo is detectable. The
production tiled sigma now consumes this generator; derivative consumers are
not migrated by this change. End-to-end speed must still be measured rather
than inferred from support generation alone.

## Generated-node tile implementation

The production cursor uses a beta-first traversal, followed by the compatible
alpha traversal for each completed beta string. This is a reordering of the
same paired-seed intersections, not a Cartesian scan followed by membership
testing. Its prefix-depth bound is $2m+2$ rather than $m+1$; the asymptotic
enumeration bound above is unchanged. A fixed-beta cursor supplies the exact
support of each exterior-transform column without scanning every alpha
index and testing it against every seed.

`DirectCiSigmaAction::apply_tiled` consumes consecutive generated nodes for
one beta string, subdividing the alpha list to fit its workspace target.
It reads that column's alpha excitation halo and only its actual beta
single-excitation neighbors. It then scatters generator channels, contracts
the integral kernel, and gathers through the reverse excitation links.
Input coefficients outside $\mathcal S$ are zero. The sigma callback may
receive output outside $\mathcal S$; the subsequent supported adjoint
transform discards those entries before projection. Thus the low-level
sigma can be compared to the full, unrestricted sigma, while the selected
structure action remains $F^\dagger H_{\mathrm{orth}}F$.

The shared numerical payload for a tile of $t$ alpha targets, excitation
halo size $h$, beta single-excitation degree $d_\beta$, and $b$ RHS is

$$
8b\,[3h+t(d_\beta+2d)]\quad\text{bytes},
$$

where $d=m(m+1)/2$ is the packed channel count. Thread workers own disjoint
rows of this shared payload; no worker has a complete beta-range reduction
array. Index maps of size $U_\alpha$ remain single shared metadata arrays.
The 64 MiB sigma target determines tile extent, not physical-space admission;
the minimum one-row tile can exceed this target for an exceptionally wide
caller block. Likewise, the panel-cache target excludes a transient incoming
panel and allows at least one panel. These targets must not be reported as a
strict process-memory cap.

The structure action expands $Pc$ once and contracts $P^T$ once per RHS batch.
Input and output panel caches expose build, projection, eviction, and retained
payload counters. If all needed panels fit, each input panel is built once,
and each H and M panel is projected once. If panels are evicted, replay is
still possible and is reported explicitly. This is not a claim of no replay
at arbitrary active-space size, nor does it remove the separate derivative
path's remaining dense-carrier limitations.

### CAS(12,12) measurements, 2026-10-15

Hanhai22 Slurm jobs 1259362 (8 cores) and 1259369 (32 cores) tested the first
1, 3, 20, and 132 singlet Rumer structures of the same synthetic integral
fixture. 132 is the complete covalent structure set. Timings include the
complete structure H/M action, not only support enumeration or sigma; they
are not molecular SCF timings. RSS below comes from native processes without
the full-carrier reference allocation. Separate parity processes compare
the entire H/M image, not merely a Rayleigh quotient.

The old executable was not rerun. Saved jobs 1259142 (8 cores, anode26) and
1259097 (32 cores, anode20) supply the baseline. New jobs ran on anode40 and
anode24, respectively; all report AMD EPYC 9654, but these are not same-node
matched timings. The anode20 request was cancelled while pending because
that node was fully occupied.

| Cores | Structures | Old cold / warm action (s) | New cold / warm action (s) | Old / new peak RSS (MiB) |
| --- | --- | --- | --- | --- |
| 8 | 1 | 0.1657 / 0.1289 | 0.2322 / 0.1976 | 52.2 / 34.8 |
| 8 | 3 | 0.2054 / 0.1578 | 0.2513 / 0.2158 | 53.8 / 35.9 |
| 8 | 20 | 0.3461 / 0.2995 | 0.2941 / 0.2523 | 58.2 / 31.8 |
| 8 | 132 | 0.1767 / 0.1369 | 0.2409 / 0.1979 | 57.9 / 32.0 |
| 32 | 1 | 0.1209 / 0.0559 | 0.3413 / 0.2695 | 79.6 / 37.1 |
| 32 | 3 | 0.1420 / 0.0642 | 0.3715 / 0.2943 | 86.4 / 33.7 |
| 32 | 20 | 0.1896 / 0.0988 | 0.4401 / 0.3481 | 88.9 / 37.2 |
| 32 | 132 | 0.1306 / 0.0602 | 0.3817 / 0.2795 | 90.1 / 37.9 |

The largest 32-core full-image errors were $1.34\times10^{-14}$ for H and
$1.91\times10^{-15}$ for M. Each RHS was expanded and projected once; each
of the four panels was built once and projected once per operator, with
zero input or output evictions. The 8-core run used 15 smaller panels with
the same once-only property. All 12 selected regression tests passed.

This implementation establishes the generated-node numerical path and reduces
measured memory, but **does not yet achieve a time improvement**. In particular,
32-core execution regresses markedly. The column-at-a-time scheduler performs
serial neighbor callbacks and several team barriers per beta column, splitting
the integral contraction into small per-worker matrix products. This is the
next scheduling issue to measure and remove, not a reason to change numerical
thresholds or assume that a smaller support automatically makes the code faster.

A bounded phase-timing run of the new binary (1259375, anode24, 32 cores)
locates 84--86% of the combined cold-plus-warm action time inside `ci_sigma`:
0.525, 0.576, 0.666, and 0.560 s for the four structure counts. These include
the scheduler's callbacks and synchronization, not just GEMM. Thus the main
remaining cost is in the generated-node sigma execution, not repeated
structure expansion/projection. The phase timer does not separate callback,
barrier, and arithmetic time; attributing their individual shares requires
finer instrumentation. The next implementation should batch independent beta
node groups into larger tasks, while retaining shared bounded storage and
exclusive output ownership.

Raw logs: `results/diagnostics/orthogonal-tiles-20261015/selected-1259362/`
and `selected-1259369/`; saved baselines: `selected-1259142/` and
`selected-1259097/` under the same directory.

### Column-batch generation and direct adjoint accumulation

The column-at-a-time path was followed by a multi-column implementation. Its
first version still generated every intermediate product serially, reconstructed
beta groups from that stream, formed a rectangular output slab, and copied the
slab into resident exterior-transform panels. Those steps are scheduling and
storage choices, not requirements of the Hamiltonian.

For either support depth, define its beta projection by

$$
\mathcal B_d=\{b:\exists a,\ (a,b)\in\mathcal S_d\},
\qquad \mathcal S_0=\mathcal S,\quad \mathcal S_1=\mathcal I.
$$

The beta prefix traversal generates exactly this projection: any surviving
seed can be paired with the earliest alpha string, whose prefix counts dominate
every seed alpha string. Complete-covalent traversal uses the corresponding
combined-prefix bound. No alpha products need to be visited to generate
the beta columns. A separate fixed-beta traversal then generates only its
compatible alpha strings. Different columns are independent and can be generated
in parallel. Only the beta list is retained, not the complete product list.

A numerical batch reserves bounded alpha slots for its beta columns and merges
their single-excitation neighborhoods. Large alpha columns are split with a
resumable cursor; a continuation does not restart the column. Each worker fills
disjoint slots of one shared index buffer. Compaction removes unused slots,
not physical directions. The planner retains the same 64 MiB slab target and
the previous conservative two-panel bound, so this change is not a larger-tile
experiment. Per-column cursor stacks contain only prefix metadata.

Let $Q=KZ$ denote the integral-contracted generator channels in a batch.
The adjoint contribution to an output coordinate is

$$
y_{ab}\mathrel{+}=\sum_{t,\lambda}
 E_{\lambda,t;ab}\,Q_{t\lambda}.
$$

The temporary scalar sum can be added directly to the caller's output panel.
Storing all such sums in an intermediate rectangular image and then adding that
image is algebraically unnecessary. The accumulator callback now performs this
addition directly; the rectangular output slab and its subsequent copy are
deleted. The output owner pins one panel until its taskgroup finishes, then may
evict it. Beta/RHS columns and disjoint alpha ranges provide exclusive writes,
without atomics or thread-private carrier arrays.

For a batch with $h$ alpha halo rows, $B$ source beta columns, $r$ intermediate
products, $d$ generator channels and $b$ RHS, numerical slab storage changes from

$$
16b(hB+rd)\quad\text{to}\quad 8b(hB+2rd)\quad\text{bytes}.
$$

Index metadata, bounded input/output panel caches, library workspaces and
allocator retention are separate; this is not a bound on process RSS.

For the selected-structure caller, $F=P_{\mathcal S}F$ implies
$F^\dagger=F^\dagger P_{\mathcal S}$. Its output panel therefore requests the
accumulator only on exact depth-zero support rows. The byte support mask is
generated once per resident panel by fixed-beta paths, shared by its RHS, and
charged to the panel-cache budget. Eviction discards both the panel and its
mask. Small cases may fit all panels in that fixed budget; this is not an
allocation sized unconditionally as $U_\alpha U_\beta$. With increasing carrier
size, residency is capped and panels are evicted. The existing minimum-one-panel
exception can exceed the target, with storage linear in $U_\alpha$ for a single
beta column. There is no separately allocated global product mask or retained
list of $\mathcal I$.
Importantly, this restriction applies to final output, not the generator-channel
intermediates on $\mathcal I$: paths leaving and returning to $\mathcal S$
remain intact. Generic sigma callers request all rows and retain the full image;
tests independently compare that image and the selected structure H/M action.

### Batch revision measurements, 2026-10-15

Slurm 1259723 passed all 12 regression tests. Jobs 1259728 (32 cores, anode20)
and 1259729 (8 cores, anode26) then ran only the new executable on the historical
benchmark nodes. Old executables were not rerun. Core placement and concurrent
node load are not matched: the historical 32-core affinity was 64--95, whereas
1259728 received 96--119,184--191. These are same-node, not identical-resource
comparisons. Native RSS is measured by GNU time; reference parity runs are
separate processes.

| Cores | Structures | Historical warm (s) | Previous batch warm (s) | Revised warm (s) | Historical / previous / revised RSS (MiB) |
| --- | --- | --- | --- | --- | --- |
| 8 | 1 | 0.12890 | 0.15363 | 0.08652 | 52.2 / 76.3 / 64.4 |
| 8 | 3 | 0.15780 | 0.16092 | 0.09286 | 53.8 / 76.7 / 65.3 |
| 8 | 20 | 0.29955 | 0.17830 | 0.11420 | 58.2 / 74.7 / 68.2 |
| 8 | 132 | 0.13694 | 0.16194 | 0.11302 | 57.9 / 71.7 / 68.0 |
| 32 | 1 | 0.05588 | 0.11438 | 0.06652 | 79.6 / 76.8 / 63.0 |
| 32 | 3 | 0.06419 | 0.14468 | 0.08477 | 86.4 / 81.8 / 72.6 |
| 32 | 20 | 0.09876 | 0.16202 | 0.09388 | 88.9 / 79.9 / 75.9 |
| 32 | 132 | 0.06015 | 0.15225 | 0.11017 | 90.1 / 82.0 / 78.1 |

The previous-batch jobs are 1259630/1259632, not the older single-column
regression 1259362/1259369. Full image absolute norm errors across the revised
8/32-core comparisons are at most 1.51e-14 (H) and 2.02e-15 (M), with unchanged
tolerances. Complete sigma images, selected structure images, unequal spin
sectors, multiple RHS, and forced panel eviction are covered by the regression
suite. No molecular convergence or energy improvement is claimed by this fixture.

The 32-core phase runs 1259633 and 1259725 (both anode15) reduce warm path/batch
preparation from 47.9 to 13.0 ms. Input reads change from 5.2 to 6.1 ms. Because
the fused writer now includes the adjoint arithmetic, compare compute+write
together: 87.7 to 81.6 ms, not the compute fields in isolation. Kernel batch
count rises from 23 to 32 and source elements read rise from 9.29 to 10.64 million;
this speedup does not come from larger tiles. Numerical slab plus target-index
peak decreases from 42,546,096 to 34,899,120 bytes (not process RSS).

The intermediate fused writer initially regressed because its column-major task
order concentrated heavy alpha-adjoint columns on too few workers. Alpha-stripe
ownership fixes that scheduling mistake. All intermediate results, including
this failed variant, are retained in `batch_comparison.csv` rather than hidden.

The revised path is faster and lower-RSS than the preceding beta-batch version,
but **does not eliminate the historical 32-core regression**: 132 structures
still take about 1.83 times the saved historical warm time. Its cold action is
0.19668 s versus historical 0.13061 s. At 8 cores all four warm actions beat the
saved historical times, but RSS remains above that older implementation.
This tradeoff and the remaining repeated beta-neighborhood reads must be kept
visible; passing parity alone does not establish a completed performance goal.

Evidence: `results/diagnostics/orthogonal-tiles-20261015/batch_comparison.csv`,
`selected-1259723/`, `selected-1259728/`, `selected-1259729/`, and
`tile-profile-1259725/`. Earlier candidate stages are identified in the CSV.

### Planned end-to-end structure-space comparison

After the selected-action implementation is validated, benchmark complete
CAS(12,12) molecular optimizations on Hanhai22, not just fixed-point H/M
actions. Use the nested structure spaces COV, ION(0,1), ION(0,2), and successive
ionic levels through FULL. Verify the generated spaces and their nesting;
do not infer structure counts from the input labels alone.

Keep geometry, AO basis, initial orbital coefficients, spin, integral backend,
optimizer, convergence criteria, thread count, and node type fixed. Retain
structure-basis Davidson throughout. Reuse existing historical results rather
than rerunning old executables; clearly label any unmatched historical setup.

For each space record:

- VB structure count and the separate alpha/beta unique-spin-string counts;
- convergence status, initial and final energies, outer iterations, HVP count,
  and Davidson action count;
- complete wall time from process start to completion and native peak RSS;
- nonoverlapping phase times for initialization, initial CI solution, orbital
  optimization, and final output, with nested CI/action timings identified as
  such rather than added twice to the total;
- source revision, job ID, node, CPU affinity, and input/output artifacts.

The primary comparison is total time to convergence and peak RSS across the
structure-space ladder. Single-action timings diagnose costs but do not replace
this end-to-end result. Different structure spaces need not converge to the same
energy. The variational minimum cannot increase on enlarging a nested space,
but a local optimizer can reach different stationary branches; record such cases
rather than treating the observed final energies as guaranteed global minima.
Unconverged or resource-limited tasks remain explicit results, not omitted rows.
This is a pending benchmark plan, not a completed performance claim.
