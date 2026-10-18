# Exact incremental admission of spin-product tiles

## Measurement boundary

The preparation timer is split into path/admission, channel layout, numerical
workspace acquisition, target addressing and beta adjacency. The path timer
includes column emission, support labels, halo discovery, budget admission and
final halo sorting; it is not a measurement of column emission alone.
Input-column cache scheduling, numerical generation, cooperative copying and
LRU insertion have separate scopes. Numerical callbacks may execute parallel
tasks even though their owner is a single coordinator. These scopes are not
all additive to the outer sigma timer. Production runs do not read clocks
at each target; the detailed scopes are opt-in diagnostics.

## Byte accounting

Let $T$ be the current target count, $h$ the alpha halo size, $b$ the number of
source beta columns and $g$ the number of target beta groups. Let $w$ denote
the RHS width, $P$ the full channel count, $s$ the channel-slab count and $p$
the number of packed slabs. Write $e=\sum_t k_t$ for the sum of the target
source-channel widths before multiplication by $w$.

Use the byte sizes $d=\operatorname{sizeof}(\mathrm{double})$,
$i=\operatorname{sizeof}(\mathrm{int})$,
$z=\operatorname{sizeof}(\mathrm{size\_t})$ and
$r=\operatorname{sizeof}(\mathrm{SigmaGatherRow})$. The existing modeled
tile payload is exactly

$$
B=aT+ce+h(ub+vg),
$$

where the action-constant coefficients are

$$
\begin{aligned}
a&=r+dw(s-p)P+\mathbf 1_{p>0}\,[i+w(i+z)],\\
c&=dwp,\qquad u=dw,\qquad v=iw.
\end{aligned}
$$

These terms account for channel slabs, halo inputs, the per-RHS image maps,
the sorted beta-gather records and the packed-layout maps. The budget does
not count referenced immutable spin graphs as per-tile numerical storage.

One new target contributes $k$ source channels and introduces
$\delta h,\delta b,\delta g$ new halo rows, source columns and target groups.
The exact change is

$$
\Delta B=a+ck+h(u\delta b+v\delta g)
 +\delta h\,[u(b+\delta b)+v(g+\delta g)].
$$

Admission therefore uses $B+\Delta B$ rather than rebuilding the complete
expression from dimensions for every target. Most targets add neither a halo
row nor a beta group/column, so only their row/channel contribution remains.
This removes repeated rectangular accounting; it does not approximate the
memory budget, increase tile capacity or change the physical space. A rejected
node is retained for the next tile, with its already-computed support metadata.
The legacy full payload evaluator and redundant cumulative channel sum are
removed. A test-only oracle recomputes all live component sizes independently
and checks equality, including packed tiles, multiple RHS and tiny budgets.

### Integer admission and saturated neighborhoods

All terms above are integer byte counts. The implementation evaluates the
coefficients and incremental candidate in unsigned 128-bit integer arithmetic,
not extended floating point. A rejected candidate can exceed the process
`size_t` budget; narrowing occurs only for the admitted payload, which is at
most that budget. With 32-bit signed width, slab, pair and per-spin dimensions,
the largest coefficient contribution is bounded by a constant times the
product of four 31-bit dimensions, hence is strictly below $2^{128}$. In the
incremental expression the stored payload is already bounded by `size_t`,
and new row/column/group counts are bounded by the existing int-indexed spin
and target spaces. This preserves exact comparison at the budget boundary
without rounded byte counts or overflowing a 64-bit rejected candidate.

Let $\mathcal H\subseteq\mathcal U_\alpha$ be the current alpha halo and
$\mathcal B\subseteq\mathcal U_\beta$ its source columns. The exact neighbor
sets of a new target can only add elements of these same single-spin universes.
Consequently,

$$
\mathcal H=\mathcal U_\alpha\Rightarrow\delta h=0,
\qquad
\mathcal B=\mathcal U_\beta\Rightarrow\delta b=0.
$$

Checking the set cardinalities is sufficient: the existing unique-index maps
already guarantee the subset relation. After saturation, querying every
neighbor again cannot change the union and is skipped. New target beta groups
still contribute $\delta g$ and are budgeted normally. This identity applies
to arbitrary selected families, regardless of whether their many-electron
space is complete. It adds no state, cache or per-worker numerical storage.

## Shared prefix certificate

For the certified complete-covalent seed family only, let $A_j$ and $B_j$ be
the alpha/beta prefix electron counts. Depth-zero membership is

$$
(a,b)\in S\quad\Longleftrightarrow\quad
A_j+B_j\ge j\ \text{for every prefix }j.
$$

The deficient-prefix set $D=\{j:A_j+B_j<j\}$ already determines the exact
source-channel envelope through its first/last elements. Thus the same
certificate gives both membership ($D$ empty) and that envelope. Packed-tile
generation previously queried membership and then scanned the same prefixes
again for frontier targets. It now computes both labels in one pass and
retains them if the node crosses a tile boundary. With packing disabled, the
ordinary membership path is unchanged. Arbitrary selected families retain
their own generator labels: their unrestricted channel envelope is **not**
used to infer membership or closure.

The changes add only action-constant byte coefficients and current-node
metadata. No product-sized cache, larger numerical tile, per-worker carrier,
new parallel region or relaxed numerical threshold is introduced. The
incremental byte identity applies to arbitrary selected structure spaces,
independently of whether the complete-covalent certificate is available.
