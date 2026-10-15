# Exact contractions for the analytic orbital Hessian diagonal

## Scope

The reference frame used to define the preconditioner is orthonormal in the
atomic-orbital (AO) overlap metric. This does not change the nonorthogonal
valence-bond (VB) orbitals, their allowed supports, the nonredundant coordinates,
or the orbital update. The results below concern construction of the analytic
fixed-structure preconditioner, not the exact Hessian-vector product (HVP).

Let $B_{\mu p}$ be the AO coefficient of reference orbital $p$, and let
$(\mu\nu|\kappa\lambda)$ denote an electron-repulsion integral. For a symmetric
AO density $D$, define

$$
J[D]_{\mu\nu}=\sum_{\kappa\lambda}
 (\mu\nu|\kappa\lambda)D_{\kappa\lambda},\qquad
K[D]_{\mu\nu}=\sum_{\kappa\lambda}
 (\mu\kappa|\nu\lambda)D_{\kappa\lambda},\qquad
G[D]=2J[D]-K[D].
$$

The consumers require only projected diagonal entries

$$
j_p[D]=\sum_{\mu\nu}B_{\mu p}J[D]_{\mu\nu}B_{\nu p},\qquad
k_p[D]=\sum_{\mu\nu}B_{\mu p}K[D]_{\mu\nu}B_{\nu p}.
$$

For inactive orbital $c$, $D^{(c)}_{\mu\nu}=B_{\mu c}B_{\nu c}$ gives
$j_p[D^{(c)}]=(pp|cc)$ and $k_p[D^{(c)}]=(pc|pc)$. These identities hold
irrespective of the selected VB structure space.

## Active sources and accepted-point reuse

For active indices $a,w,z$, let $\gamma$ and $\Gamma$ denote the accepted
one-electron density and the integral-symmetric two-electron density
representative, respectively, in the reference active frame. Define symmetric
source weights

$$
X^{(a)}_{wz}=\frac12\left(\delta_{wa}\gamma_{az}
                           +\delta_{za}\gamma_{aw}\right),
$$

$$
Y^{(a)}_{wz}=\frac12\left(
 \Gamma_{awza}+\Gamma_{waza}+\Gamma_{azwa}+\Gamma_{zawa}\right).
$$

The corresponding AO density is $D[W]=B_{\rm act}WB_{\rm act}^{T}$.
The total one-electron source satisfies the exact identity

$$
\sum_a X^{(a)}=\frac12(\gamma+\gamma^T).
$$

Consequently its operator image is the sum of the $X^{(a)}$ images; it does
not require an additional AO contraction. An untruncated Householder QR of
the packed $2n_{\rm active}$ sources gives at most
$\min[n_{\rm active}(n_{\rm active}+1)/2,2n_{\rm active}]$ basis densities.
This is an exact change of source basis, not a numerical rank truncation.

The accepted transformed integrals already contain $\mathcal J_{p,wz}=(pp|wz)$.
Thus, for every active source $W$,

$$
j_p[D[W]]=\sum_{wz}\mathcal J_{p,wz}W_{wz},\qquad
g_p[D[W]]=2j_p[D[W]]-k_p[D[W]].
$$

In the output-owned row layout, only exchange needs a new active-source
AO-integral contraction. Coulomb is reconstructed in the small active-pair
space. The original mean-field consumers are preserved
by this exact reconstruction; no structure-response term is added to the
analytic diagonal.

## Exclusive output-row contraction

Assign AO output row $\mu$ to one worker. For each ordered $j$, traverse the
existing symmetric compressed-row (CSR) integral row for pair $(\mu,j)$.
For its right pair $k\ge l$ with raw integral value $v$,

$$
J[D]_{\mu j}\mathrel{+}=v(2-\delta_{kl})D_{kl},
$$

$$
K[D]_{\mu k}\mathrel{+}=vD_{jl},\qquad
K[D]_{\mu l}\mathrel{+}=vD_{jk}\quad(k\ne l).
$$

No unique-integral half weights belong in this directed-CSR contraction.
Coincident AO indices follow the same expressions. Workers never update the
transposed row, so neither atomic scatters nor complete thread-local AO
responses are needed. Density channels are contiguous for SIMD evaluation.

Density symmetry is also a locality identity, not only a storage symmetry:

$$
D_{jl}=D_{lj},\qquad D_{jk}=D_{kj}.
$$

With $j$ fixed, the exchange inputs are read from one contiguous AO density
column. A strided access through successive density columns is unnecessary.
The Coulomb orientation likewise follows the packed right-pair ordering.
All channels use the same shared density data; no private density copy is made.

Inactive Coulomb rows are projected before being discarded:

$$
j_p[D^{(c)}]=\sum_\mu B_{\mu p}
 \left[\sum_j B_{jp}J[D^{(c)}]_{\mu j}\right].
$$

Each worker accumulates a small orbital-by-channel panel. This removes the
separate core AO-pair integral action and its AO-pair-by-core result. Exchange
uses one shared AO response batch; one shared contiguous AO image and disjoint
MO-column product panels give its projected diagonals.

## Work and memory

Let $N$ be the AO dimension, $T$ the worker count, $b$ the density batch width,
$r$ the MO projection-panel width, and $E$ the number of unique stored ERIs.
The incremental contraction workspace scales as

$$
O(N^2b+N^2+TNb+TNr),
$$

instead of $O(TN^2b)$. Batch width and worker count are admitted together
under the same workspace budget, including sources, response, AO staging,
projection panels, density-building temporaries and small outputs. The large
integral graph and accepted tensors are shared and are not part of this
incremental workspace budget. This budget accounts for explicit application
buffers; it is not a strict heap or process-RSS cap. Allocator overhead and
internal Eigen/BLAS packing workspaces are not included. For exceptionally
large AO dimensions, the unavoidable single-channel/single-worker working
set is the lower bound, rather than enabling additional large replicas.

An ordered-row batch reads approximately four times as many integral records
as a canonical unique-integral batch, owing to pair exchange and left-pair
orientation. Its work remains $O(EC)$ for $C$ density channels; this is not
an asymptotic reduction of exact four-index integral contraction. Wider batches
can reduce integral/metadata rereads, while accepted Coulomb reuse and fused
core projection eliminate separate calculations. Both benefits must be measured
against the extra CSR orientation reads and projection work.

For admission, let $L$ be the number of directed CSR entries, and let $L_d$
count entries in rows whose left AO pair has equal indices. An ordered-row
batch reads $2L-L_d$ entries. If the admitted row and canonical batch widths
are $b_r$ and $b_s$, respectively, the main-contraction record counts are

$$
R_r=(2L-L_d)\left\lceil\frac{C}{b_r}\right\rceil,\qquad
R_s=E\left\lceil\frac{C}{b_s}\right\rceil.
$$

Let $L_2$ be the runtime-reported per-core L2 cache capacity. The packed
canonical output footprint per worker is

$$
P_s=8\frac{N(N+1)}{2}b_s\quad\text{bytes}.
$$

The row layout is admitted when $R_r<R_s$ and $P_s>L_2$, provided $L_2$ is
available. Otherwise the canonical layout is retained. This combines an
integral-record traffic estimate with a cache-capacity heuristic; it is not
a complete wall-time model or a guarantee of optimal layout. The cache
test accounts for the packed output, not the entire density/integral working
set, and reported cache topology can vary across platforms. The saved core
action is deliberately not priced as a free benefit because its fused
projection also costs work. No molecule-specific switch or change of
numerical operator is involved. Small operators can retain efficient
canonical packed scatter even when a wider row batch reduces record reads.

The canonical integral stream remains a required representation when the CSR
edge address space is exceeded. The row-owned method does not construct a
second integral graph or silently substitute an approximate operator.
