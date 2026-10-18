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
