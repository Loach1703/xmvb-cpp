# Exact channel support in a selected-structure direct-CI action

## Source support is not channel support

Let $\mathscr A$ map selected VB structure coefficients into the orthogonal
alpha/beta spin-string carrier. Its application is generated in tiles; the
matrix is not assembled. For any selected coefficient vector $v$,

$$
X(v)=\mathscr A v,\qquad \operatorname{supp}X(v)\subseteq\mathcal S.
$$

Here $\mathcal S$ is a certified symbolic support, not a numerical amplitude
cutoff. In the current upper-triangular exterior representation it is generated
from the selected source products. It may conservatively include entries that
cancel for a particular wave function.

For packed orbital pairs $P=(p,q)$ with $p\geq q$, define the spin-summed
one-body generators by

$$
E_{pp}=n_{p\alpha}+n_{p\beta},\qquad
E_{pq}=\sum_{\sigma=\alpha,\beta}
\left(a^\dagger_{p\sigma}a_{q\sigma}
+a^\dagger_{q\sigma}a_{p\sigma}\right),\quad p>q.
$$

For an intermediate spin-string pair $r$, the channel vector is

$$
F_r(v)_P=(E_P X(v))_r.
$$

The existing Hamiltonian kernel $K$ includes the absorbed one-electron and
normal-ordering contributions, so the selected bilinear form is

$$
u^T H_{\mathrm{VB}}v=\sum_r F_r(u)^T K F_r(v).
$$

Only the one-generator closure of $\mathcal S$ can contribute to this sum.
Within that closure, however, the channel dimension need not be the complete
$d=n_{\mathrm{act}}(n_{\mathrm{act}}+1)/2$.

## A structural zero certificate

Define

$$
\mathcal J_r=\left\{P:\exists t\in\mathcal S,
\ (E_P)_{rt}\ne0\right\},\qquad k_r=|\mathcal J_r|.
$$

Then $F_r(v)_P=0$ for every $v$ whenever $P\notin\mathcal J_r$. This follows
from the operator's matrix-element support, not from current coefficient sizes.
For a diagonal channel, $r$ itself must be in $\mathcal S$ and site $p$ must
be occupied in at least one spin sector. For an off-diagonal channel, at least
one allowed single-spin replacement must have its preimage in $\mathcal S$.
No overlap inverse, cofactor or determinant-pair table is needed for this test.

## Universal image versus selected bilinear form

To preserve a **full-image** sigma contract, retain every output channel:

$$
W_r=K_{:,\mathcal J_r}F_r(v)_{\mathcal J_r}.
$$

Computing only $K_{\mathcal J_r,\mathcal J_r}F_r$ is not a valid replacement
for that full-image contract. It is sufficient only when the complete left
consumer is known to have the corresponding certified support. In particular,
the exact selected bilinear form may be written as

$$
F_r(u)_{\mathcal J_r}^T
K_{\mathcal J_r,\mathcal J_r}
F_r(v)_{\mathcal J_r}.
$$

This identity assumes neither completeness nor invariance of the selected
VB structure subspace: both left and right maps remain $\mathscr A$.

For the integral cotangent of this bilinear form,

$$
G=\sum_r F_r(u)F_r(v)^T,
$$

each row's contribution is confined to $\mathcal J_r\times\mathcal J_r$.
Exchange-orbit multiplicities remain exactly those of the existing quotient.

## Derivatives and scope

The symbolic upper-triangular support is fixed by the selected source family,
not by zeros of the accepted orbital transform. Derivatives of its exterior
polynomials have the same symbolic support. Thus orbital tangents, CI tangents
and higher jets can use this certificate as long as their complete source maps
use that same family. A derivative must **not** reuse a mask obtained by
discarding currently zero but structurally allowed amplitudes.

Skipping a multiply-add whose supplied scalar is exactly zero is separately
algebraically valid for a forward multiplication. It is not a certificate for
deleting a derivative edge. Current production derivative formulas remain the
reference; the census/prototype does not change them.

For $R$ intermediate nodes, dense channel products cost $O(Rd^2)$. A universal
full-output sparse product has arithmetic count $O(d\sum_r k_r)$, while the
integral cotangent has count $O(\sum_r k_r^2)$. These are sparse-work bounds,
not guaranteed wall-time gains or a new power law in active-space size. They
do not remove the exponential carrier growth.

## Which orbital cotangents can use a projected image?

For the current upper-triangular orbital factor $R$, its relative direction

$$
U=\delta R\,R^{-1}
$$

is upper triangular. The corresponding exterior tangent is a one-body
generator action $\widehat UY$, plus a transformed selected CI tangent.
The symbolic support is downward closed under these upper generators, so
both tangents remain in $\mathcal S$. This is differentiation of the
nonorthogonal representation, not an orbital-rotation update.

Let $D$ be a generator cotangent with
$\delta\mathcal L=\langle D,\delta R R^{-1}\rangle_F$. Its coefficient
pullback, restricted to upper directions, is

$$
g_R=\operatorname{upper}(D R^{-T}).
$$

Because $R^{-T}$ is lower triangular, an upper output entry with $i\leq j$
uses only $D_{ik}$ with $k\geq j\geq i$. Strictly lower entries of $D$ cannot
contribute. The directional pullback has the same property:

$$
\delta g_R=\operatorname{upper}
\left(\delta D R^{-T}-D R^{-T}\delta R^T R^{-T}\right).
$$

All right multipliers are lower triangular. This supports using projected
Hamiltonian images for these upper-generator cotangents: their right tangent
vectors already lie in $\mathcal S$.

The current streamed generator implementation explicitly skips $p>q$ when
`symmetric_sum=false`; product-jet requests use that setting and untransposed
right arguments. By contrast, a request for `upper(D+D^T)` or a transposed
generator cotangent can require information from lower generators. It must
not inherit the shortcut merely because the final returned matrix is upper.
The existing tests already distinguish that orientation issue. The new
projected operator contract must name its consumers, not silently redefine
the general full-image sigma interface.

## Direct boundary-channel generation for the complete covalent family

Number sites from 1 to $n$, and define the combined prefix deficit of a node by

$$
d_\ell=\ell-\sum_{i=1}^{\ell}(n_{i\alpha}+n_{i\beta}).
$$

For the complete covalent family's upper-triangular source support,
$d_\ell\leq0$. Its one-generator closure satisfies $d_\ell\leq1$. At a
boundary node let $f$ and $l$ be the first and last positions where $d_\ell=1$.
An allowed preimage must move one electron of a fixed spin from an occupied
site $q$ to an empty site $p$. This increases all prefix counts on $[p,q)$
by one and leaves every other prefix unchanged. It repairs every deficit
**if and only if**

$$
p\leq f,\qquad q>l.
$$

Necessity follows by covering both extreme deficient prefixes; sufficiency
follows because all deficiencies are at most one and all other prefixes can
only improve. Define the early-empty and late-occupied sets

$$
\mathcal E_\sigma=\{p\leq f:n_{p\sigma}=0\},\qquad
\mathcal O_\sigma=\{q>l:n_{q\sigma}=1\}.
$$

The complete boundary channel set is therefore

$$
\mathcal J_r=\operatorname{pack}
\left[(\mathcal E_\alpha\times\mathcal O_\alpha)
\cup(\mathcal E_\beta\times\mathcal O_\beta)\right].
$$

There are no diagonal channels at a boundary node. If $e_\sigma$ and
$o_\sigma$ denote these set sizes, the channel count has the closed form

$$
k_r=e_\alpha o_\alpha+e_\beta o_\beta
-|\mathcal E_\alpha\cap\mathcal E_\beta|
 |\mathcal O_\alpha\cap\mathcal O_\beta|.
$$

Find the deficient interval in $O(n)$ work, then enumerate the two bit-set
Cartesian products. They emit at most $2k_r$ entries; duplicates can be removed
while emitting the beta contribution. Channel-index generation consequently
costs $O(n+k_r)$, without a separate source-membership query for every orbital
pair or spin link. A dense bitmap's initialization would have its own cost;
an emitted-index interface need not allocate or scan all $d$ entries.

The diagnostic uses two bitmap words only at n=4 and n=12, and compares the
emitted set against the general preimage definition. This certificate applies
to the **complete covalent source family**, including unequal fixed spin
populations with $n_\alpha+n_\beta=n$. It must not be applied to an arbitrary
selected or ion-truncated family merely because that family is small. Such
families still use their actual generated source/preimage support. Both cases
belong to the same exact $\mathcal J_r$ contraction, not different physical
approximations or memory-triggered model choices.

Hanhai22 job 1265005 verifies the bit-set construction on all 191268 boundary
nodes of the n=12 covalent census, with no mismatch against the general
preimage-based channel mask. It also checks n=4 unequal-spin boundary cases
against production generator basis actions. The final regression suite passes
11/11. This validates the channel certificate, not a production projected
sigma or a wall-time speedup.

## Hybrid interior/boundary product

An explicit source-projected image contract permits a hybrid contraction.
Partition each generated tile into interior nodes $r\in\mathcal S$ and
boundary nodes $r\in\mathcal I\setminus\mathcal S$, without moving or
duplicating any carrier. Interior channel rows retain the contiguous dense
product. Boundary rows evaluate only

$$
W_{r,\mathcal J_r}
=K_{\mathcal J_r,\mathcal J_r}F_{r,\mathcal J_r},\qquad
W_{r,\overline{\mathcal J_r}}=0.
$$

The second equality specifies a representative for the projected consumer,
not the full Hamiltonian image. It is valid because that consumer's left
channel vector vanishes on $\overline{\mathcal J_r}$. For a directional
kernel, use the same channel sets in

$$
\dot W_{r,\mathcal J_r}
=K_{\mathcal J_r,\mathcal J_r}\dot F_{r,\mathcal J_r}
+\dot K_{\mathcal J_r,\mathcal J_r}F_{r,\mathcal J_r}.
$$

Thus the channel-product work per right-hand side becomes

$$
O\!\left(|\mathcal S|d^2+
\sum_{r\in\mathcal I\setminus\mathcal S}k_r^2\right),
$$

in addition to support generation, scatter and gather. This hybrid bound
leaves interior structural zeros unexploited. Leaving the integral Gram
contraction dense also retains its $O(|\mathcal I|d^2)$ work; accelerating
one product is not a claim that the whole derivative has this reduced cost.

The boundary row requires at most $d$ channel indices and $d$ copied input
scalars, permitting in-place products without copying the channel slab.
Parallel row ownership therefore requires $O(Td)$ small scratch for $T$
workers, not $T$ copies of the spin-string product carrier. The source/halo
partition is mathematical and independent of the numerical tile budget.

For a response/state pair $F_r(v),F_r(u)$, the same emitted channel set
also gives the exact boundary contribution

$$
G_{\mathcal J_r,\mathcal J_r}\mathrel{+}=
w_r F_r(v)_{\mathcal J_r}F_r(u)_{\mathcal J_r}^{T},
$$

where $w_r$ is the exchange-orbit multiplicity. This contribution must be
consumed **before** an in-place kernel product overwrites either argument.
Each row owner may therefore add its Gram contribution and form its required
images in the same boundary visit. Worker-local Gram matrices already used
by the dense reduction need not be duplicated. Interior rows retain batched
matrix products, and the existing Gram reduction completes before its scratch
is reused for an in-place interior image product. Neither a physical response
nor a channel with a permitted preimage is removed by this fusion.

This is an exact contraction identity, not a claim of measured acceleration.
The C12 projected-kernel implementations did not beat the saved native
control and were removed from production. Their numerical and timing evidence
is recorded in `results/diagnostics/c12-projected-channels-20261017/`.

## C12 census before implementation

Hanhai22 job 1264972, CAS(12,12), even-exchange representatives, $d=78$.
No production SCF run or old native baseline was repeated. The single-structure
source is one product of six disjoint singlet pairs (64 original products),
not the complete covalent family.

| Symbolic family | Intermediate nodes | Mean $k_r$ | $\sum k_r/(Rd)$ | $\sum k_r^2/(Rd^2)$ |
|---|---:|---:|---:|---:|
| One covalent structure | 261868 | 18.5752 | 0.23814 | 0.11333 |
| Complete covalent family | 304590 | 24.8429 | 0.31850 | 0.17335 |
| Full source carrier | 427350 | 61.3427 | 0.78644 | 0.61951 |

At 32768 diagnostic rows per tile, equal-mask groups average only 3.37, 1.99
and 1.17 rows, respectively. A naive equal-mask GEMM-per-group is therefore
not justified by the mean channel savings alone. A large group does occur in
the sparse families; the complete grouping distribution must be retained.

Further census 1264993 separates source nodes from the boundary of their
one-generator closure:

| Family | Source nodes / mean channels | Boundary nodes / mean channels |
|---|---:|---:|
| One covalent structure | 68253 / 48.8445 | 193615 / 7.9047 |
| Complete covalent family | 113322 / 50.9439 | 191268 / 9.3787 |
| Full source carrier | 427350 / 61.3427 | 0 |

The covalent family's boundary is 62.8% of its intermediate nodes but requires
only about 9.38 channels per node. This supports investigating an exact
interior/boundary contraction, using the existing source-support label rather
than a molecule-specific threshold or numerical-amplitude screening.

There is an important interface gate before such an implementation: a public
full-image sigma must continue to return its complete image. A projected
boundary kernel can restrict both channel axes only after every left consumer
has supplied its support/cotangent contract. Generator moments, their triangular
pullbacks, CI forcing and adjoints must all be audited; projection cannot be
inferred from a memory limit. A selected congruence need not be closed under
the ambient Hamiltonian, but its left and right maps must remain exact.

The first numerical probe (1264978) retains full output and skips only exactly
zero inputs. It agrees with GEMM to about 2e-17 relative error, but is slower:
1.34x for one structure, 1.64x for COV and 2.94x for FULL. Its row-wise SAXPY
loop forfeits GEMM's reuse; fewer nominal operations alone are not a speedup.
This prototype is archived, not installed in a production path. Register/block
reuse or an explicit projected-consumer contraction still needs an efficient
implementation and an end-to-end qualification.

The census streams the support cursor and retains masks only within one tile.
Its global histograms have fixed size; there is no $U_\alpha U_\beta$ mask
allocation. A small n=4 complete-carrier oracle, including unequal spins and
an empty spin sector, compares the certificate against independent production
one-body generator basis actions. That oracle is not a large-space algorithm.
