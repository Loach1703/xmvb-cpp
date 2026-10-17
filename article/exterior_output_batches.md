# Exact exterior identities and bounded output batching

## Identity factors

At a fixed orbital point an elementary exterior shear is

$$
F_j=I+c_jN_j,\qquad N_j^2=0.
$$

Here $c_j$ is a scalar orbital-elimination coefficient and $N_j$ replaces
one occupied spin orbital by an unoccupied one with its fermionic sign.
Its source/target pairs are disjoint. The complete exterior transform is
an ordered product of these shears and a diagonal occupation scale; no dense
compound matrix is required.

If $c_j=0$ **exactly**, $F_j=I$. Its primal forward, adjoint and inverse
operations may therefore be omitted. No tolerance is used. This does not
justify deleting the edge from the transform topology, because

$$
\delta F_j=(\delta c_j)N_j
$$

can be nonzero at $c_j=0$. Directional transforms keep all original edges and
their directional coefficients. Regression includes a finite difference at
the identity transform with a nonzero orbital direction.

The retained supported transform keeps the existing column-local numerical
traversal, exact support mask, ordered shears and per-worker byte mask. It
creates no new numerical workspace. Sharing support predicates across RHSs
was tested but did not establish an end-to-end advantage over this simpler
column-local loop; that alternative is not retained.

## Budgeted output task groups

An output writer adds independent beta/RHS/alpha-stripe contributions to
caller-owned compact panels. One task group per panel is not mathematically
required. Let $m_p$ charge the numerical entries and CSR metadata of panel
$p$, and let $M$ be the existing output-panel budget. Batch consecutive panels
subject to

$$
\sum_{p\in\mathcal B}m_p\leq M.
$$

Borrow the panels serially, then run the existing disjoint updates across the
whole batch. A normal LRU lookup makes each borrowed panel newer than every
unborrowed resident panel. Since the complete borrowed batch fits, adding the
next member can evict only unborrowed older entries. List splicing does not
invalidate references. The task group completes before the next batch mutates
the cache. A single unavoidable oversized panel follows the cache's original
one-panel lifetime rule, never coexisting with a borrowed second panel.

Batching therefore needs only pointers and column indices, not copied
numerical panels, a larger cache budget or more workers. An evicted partial
image is consumed through the same linear projection as before. Arbitrary
selected VB structure spaces retain their original exact projection; no
selected-space closure is assumed.

This removes repeated task-group scheduling/completion and exact identity
operations, not Hamiltonian terms, CI responses, projections or precision
requirements. It is not a lower asymptotic physical-HVP scaling claim.
Tiny-budget eviction tests accompany the resident case; a resident-only
benchmark cannot establish safe panel lifetimes. Native complete-convergence
measurements, not the algebra alone, determine the end-to-end benefit.
