# Two-channel-slab directional sigma: prototype assessment

Status: the single-direction prototype passed independent numerical tests but
was not retained. Hanhai22 jobs 1264662/1264664 took 31.462584/31.070410 s
end-to-end, versus the saved 30.863387 s baseline. Both retained 4 steps and
36 physical HVPs. The single-direction production path therefore keeps ordinary
GEMM and the [lifetime-exact halo bound](sigma_halo_budget.md). The temporal
CI-then-orbital schedule described below has since been implemented and qualified
through the physical coupled action; see [temporal sharing](temporal_ci_responses.md).
Its benefit comes from eliminating a second accepted scatter and graph sweep,
not from assuming that an in-place single-direction GEMM is intrinsically faster.

## Required operations

Let $\mathcal E$ be the one-body channel scatter, $Y$ the accepted carrier and
$Z$ its orbital/coefficient direction. The channel matrices are
$A=\mathcal E Z$ and $B=\mathcal E Y$. Let $K$ be the Hamiltonian channel
kernel, $\dot K$ its integral direction, and $W$ the diagonal matrix of exact
spin-exchange orbit multiplicities. Before overwriting either input, form the
integral cotangent $G=A^T W B$. The directional channel images are

$$
U=A K^T,\qquad V=B\dot K^T,\qquad C=U+V.
$$

The previous kernel kept A, B and C simultaneously, then reused A for V.
However each output row of either GEMM depends only on that input row. After
both scatters have completed, the coefficient halo is dead until the next
tile read. It can serve as shared micro-panel scratch:

1. Compute a bounded row block of A K^T into a disjoint slice of the dead halo.
2. Copy that block back into the same rows of A; other A rows are independent.
3. Restore B's original, unweighted exchange convention, then similarly
   compute B dot(K)^T into that halo and copy back into B.
4. Add B into A. A is the combined image, B the optional forcing image.
5. Consume both transpose gathers before the next halo read.

For the usual request that consumes only C, materializing V is unnecessary.
For each owned row block, accumulate both GEMMs directly into scratch,

$$
Q\leftarrow A_{\mathrm{block}}K^T,\qquad
Q\leftarrow Q+B_{\mathrm{block}}\dot K^T,\qquad
A_{\mathrm{block}}\leftarrow Q.
$$

This removes the separate B-to-V copyback and final A+=B traversal. If the
caller actually requests V as an independent output, both overwritten slabs
remain necessary. Adjoint-only requests need neither Hamiltonian GEMM. These
are output-dependency decisions, not numerical or cost thresholds.

No full third channel image or numerical transpose is needed. All permanent
channel views remain row-major, retaining the established GEMM/gather layout.
Scratch is shared and partitioned by worker; no worker allocates a carrier.
Each partition must hold at least one channel row. Tiny halos therefore need
a shared minimum of active_workers*d doubles, not a divide-by-zero shortcut.
The planner must explicitly bound max(halo_size, active_workers*d)+2*channel_size
and the target-map metadata. The existing 64 MiB cap must not change.

Workers own disjoint rows, so their scratch->A/B copies cannot overwrite another
worker's live input row. Existing scatter/Gram and final-writer barriers must
remain. Scratch row count follows from the dead-halo capacity and actual team,
not a tuned molecule-specific block size. Additional micro-GEMMs and copyback
costs can offset larger-batch savings; only converged native timing qualifies it.

## Why this precedes combined orbital/CI actions

At a coupled point, the same accepted B is scattered separately for orbital
Z_p and CI Z_q. Sharing B is not free if simply done with a third input halo
or more channel slabs: that can fragment batches and erase read savings.

Use temporal order instead. Within one shared generated-node tile, scatter B
once, process the CI direction first, consume its image, then reuse its A slab
for the orbital direction. The orbital direction is last, so B can become its
forcing image without harming a subsequent RHS. Keep the two integral
cotangents and two projected coupled images separate: the coarse preconditioner
uses their separate orbital/CI images, so summing them into one direction would
change the optimizer and is not an acceptable contraction-only optimization.

Nonresident points also require HY. Its channel image can be consumed before
the two direction stages, while the appropriate directional readers supply
the two reverse-generator arguments. Their input/output/projection LRUs must
share the existing caller budget, rather than duplicating a 64 MiB cache for
each direction. Arbitrary selected structure spaces remain exact; no closure
or small-residual assumption is available.

## Qualification sequence (historical design)

First implement and independently verify the single-direction two-slab kernel,
including widths>1, forcing on/off, adjoint-only requests, unequal spin spaces,
both exchange parities, small halos and nested teams. Test only the new code
on Hanhai22 with CAS12/LOFLEA, recording energy, steps, HVP, end-to-end wall and
RSS against saved results. Do not rerun old binaries. Only then generalize the
common kernel to temporal multi-RHS processing and connect both physical models;
do not create a second separately maintained contraction implementation.
