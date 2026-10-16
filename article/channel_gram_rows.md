# Row-partitioned channel Gram contraction

## Necessary contraction

For a generated spin-string tile, let $A=\mathcal E Z$ be the response
channels and $B=W\mathcal E Y$ the accepted channels, with the exact certified
exchange multiplicity $W$ already applied. The integral cotangent requires

$$
G=A^T B.
$$

There is no assumption that $A=B$ or that $G$ is symmetric. Both independent
orbital and CI cotangents remain necessary. In particular, neither the Gram
product nor its out-of-selected-space one-body halo may be dropped.

Partition the common row index into disjoint intervals $I_t$. Then

$$
G=\sum_{t=0}^{P-1} A_{I_t,:}^T B_{I_t,:}.
$$

This is an exact redistribution of the summation, including the exchange
weights and every RHS. It replaces a grid of small output products, each of
which scans every input row, with row-local products and a small output
reduction. Floating-point summation order changes; numerical equivalence must
be checked at the unchanged accuracy.

## Bounded storage from lifetimes

The coefficient input halo is dead after both scatters have finished. The
channel slabs $A$ and $B$ are separate and remain live. Reuse the dead halo
for disjoint partial Gram tiles, not for copies of either carrier.

Let $s$ be the number of available doubles in that halo/scratch region, $r$
the channel row count, $d$ the channel width and $T$ the active team size.
The row-partition geometry is used when full local Grams fit:

$$
P=\min(T,r)>1,\qquad Pd^2\leq s.
$$

Each worker borrows a full $d\times d$ partial, and one row-local product per
worker is sufficient. Otherwise the same contraction assigns disjoint output
blocks directly to workers, with no partial storage or cross-worker reduction.
For one worker this is just the necessary GEMM. This is workspace geometry,
not a gate on the physical representation or precision. Even a one-value halo
remains valid and does not create a barrier for every scalar output entry.

For a row partition: finish all partial products, synchronize, let exclusive
output-entry owners reduce the partials, then synchronize before reusing the
scratch. The final synchronization also precedes directional GEMMs and the
next input-halo refill. The common accepted channels are not overwritten by
this contraction. No new allocation, larger budget or thread-private carrier
is needed; even at large $d$, the partials stay inside the existing64MiB
double workspace.

## Work and data access

The arithmetic remains $2rd^2+O(Pd^2)$ floating-point operations. When the
full Gram fits, a square $q\times q$ output grid logically reads the channel
operands $q$ times, whereas the row partition reads each row in one local
matrix product and writes/reduces $Pd^2$ partial values. The nominal operand
volume changes from $O(rdq)$ to $O(rd+Pd^2)$. These are logical matrix-product
operands, not measured DRAM traffic: cache reuse and BLAS packing affect the
actual hardware traffic. This is not a lower polynomial order in active-space
size, nor a guarantee of end-to-end speedup.

All Hamiltonian images, forcing gathers, integral pullbacks and structure
projections are unchanged. Qualification uses the existing independent
direct-CI image/cotangent tests, tiny-halo and exchange cases, then native
C12/LOFLEA convergence on Hanhai22 against saved old measurements only.

## Qualification

The final implementation in Hanhai22 job1264786 passes all nine regression
targets. Added eight-orbital, three-electron-per-spin tests explicitly request
two and four workers with two distinct response RHS blocks. Their halo fits
full36-square partial Grams; images and integral cotangents agree with
independent non-tiled operators at the unchanged3e-11 test tolerance for
exchange parities0,+1,-1. Tests also retain strict selected support, unequal
spin spaces, temporal CI/orbital requests, zero directions, forcing-only and
adjoint-only requests, nested teams and explicit four-worker minimal halos.

The complete native CAS(12,12) LOFLEA task takes24.205419s versus the saved
27.623396s predecessor, with the same four accepted steps and36 HVPs. Final
energy is $-399.106646267332$ hartree and GNU-time peak RSS2650.949MiB.
In the separate phase run the Gram scope falls2.661783 to0.563851s and the
complete coupled-image scope9.676411 to7.989973s. C12 uses32 row partitions
of78-square Grams, borrowing1,557,504 bytes from the already allocated halo.
The number of sweeps, generated batches and source read elements is unchanged.
No old executable was rerun; raw evidence is under
`results/diagnostics/c12-row-gram-20261017/`.
