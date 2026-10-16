# Temporal sharing of accepted-state CI channels

Status: implemented through the physical orbital HVP, CI tangent operator and
CIAH coupled action. Hanhai22 job 1264713 completed the exact-integral CAS(12,12)
LOFLEA regression in four accepted steps and 36 orbital HVPs. Single and
multiple response requests use the same generated-graph contraction engine.

## Independent response outputs

At one accepted point let $Y$ be the accepted spin-string carrier, $Z_j$ an
orbital or CI direction and $\mathcal E$ the one-body channel scatter. Define

$$
B=\mathcal E Y,\qquad A_j=\mathcal E Z_j,\qquad
G_j=A_j^T W B,\qquad
C_j=A_j K^T+B\dot K_j^T.
$$

$K$ is the Hamiltonian channel kernel and $W$ contains the exact orbit
multiplicities of a certified alpha/beta exchange parity. Without an exchange
reduction $W=I$. A CI-only request has $\dot K_j=0$. Each $G_j$ independently
pulls back to that direction's one- and two-electron integral cotangents; each
$C_j$ independently gathers to its image. They must not be summed: CIAH's
coarse preconditioner retains separate orbital and structure images.

For a generated tile, scatter the accepted state once and retain $\widetilde
B=WB$. For each response, reuse the same response slab and output slab in
sequence. Consume that response's gather before overwriting either slab.
The topology, halo address maps and accepted scatter are independent of $j$.
No full spin-product coefficient table or per-worker carrier is required.

For an integral direction followed by another request, the accepted slab must survive.
Since W acts on rows and $\dot K_j$ acts on columns, form the forcing as

$$
B\dot K_j^T=W^{-1}\bigl(\widetilde B\dot K_j^T\bigr).
$$

On a representative exchange orbit, $W^{-1}$ is exactly 1 or 1/2. Apply it
to the temporary forcing image, not to the retained accepted channels. This
keeps later Gram contractions correct, including a later adjoint-only request.
For a single request the established unweighted-state GEMM order is unchanged.
Requests with no image or forcing consumer omit the unused Hamiltonian GEMM.

## Ordered responses and storage

The physical coupled action needs a CI request followed by an orbital request.
This order gives a stronger lifetime identity. After the CI Gram contraction,
overwrite its response channels with $A_qK^T$, consume the image, and overwrite
the same slab with $A_p$. The common $WB$ remains intact. After the final Gram,
no later request needs $B$: unweight it and overwrite it with $B\dot K_p^T$.
Add this forcing to $A_pK^T$ and consume the two images separately.

The row-local channel products use the already consumed input halo as scratch.
Each worker owns a disjoint range; neither a carrier nor a channel slab is
replicated per thread. This two-slab schedule is used only when every request
before the last is CI-only. An intermediate integral direction still requires
three slabs to preserve the accepted channels for later Gram products.

For $r$ generated targets, $h$ halo rows, $b$ halo columns, $g$ beta groups,
$w$ RHSs, $d$ orbital-pair channels and $t=\min(N_{\mathrm{threads}},rw)$,
the ordered schedule plans the following live payload:

$$
8\left[\max(whb,td)+2wrd\right]+4hg
\leq64\times2^{20}\ \mathrm{bytes}.
$$

The general three-slab schedule instead uses $8w(hb+3rd)+4hg$ bytes.
These bounds cover the shared double workspace and target map, not every
container in the process. Spin-link topology, small address arrays and the
caller's exterior-transform panel caches are separate. The existing caller
cache instances and their individual budgets are shared across requests,
rather than duplicated. There are still separate small $d\times d$
cotangents, so total storage is not independent of that number. Neither input
nor output callbacks may retain the ephemeral shared tile views.

The Gram and Hamiltonian products remain necessary for every independent
direction. Only accepted reads/scatters, topology construction and their
synchronization are eliminated. This changes the number of sweeps, not the
Hamiltonian's formal polynomial order. A three-slab prototype lost the common
read savings to smaller batches; it is not the schedule used for the physical
CI-then-orbital pair. Isolated kernel timing is variable and is not a substitute
for the full convergence measurement.

## Physical integration

The CIAH adapter prepares raw horizontal CI coefficients before the orbital
integral-direction sweep. A common product-jet session supplies response readers
and separate output consumers for each selected state. In an accepted-image
cache miss, the packed inputs are $[Y,\dot T_p C]$ and $[0,Tq]$ with the common
right input $[0,Y]$. The accepted Hamiltonian image is generated only by the
orbital request. On a cache hit, both requests have width one. There is no
assumption that an arbitrary selected VB space is closed under the exterior
transform $T$.

The CI output supplies $Hq$, $Mq$, integral cotangents and exterior-generator
cotangents. The existing tangent finisher applies the same state projection,
energy shifts and normalization terms without invoking sigma again. A pure
CI direction skips the orbital integral transform. An exactly zero CI direction
does not launch an extra product-jet request. The existing spin-factor plan,
when algebraically cheaper, remains a distinct exact contraction representation.

The downstream physical operations remain unchanged: per-state energy shifts,
weights, left exterior-transform derivatives, state gauge multipliers, selected
space projection, orbital adjoint pullback and the two coarse-model images.
No relaxation of precision, response equation or outer stopping rule is involved.

## Qualification

Unit tests compare the shared contractions against independent separate calls
for exchange parities $0,+1,-1$, unequal spin spaces, nonclosed selected spaces,
multiple states, zero weights and resident/nonresident accepted images. The
CAS(12,12) operator comparison has relative differences below $10^{-15}$.
The physical full-chain test also compares mixed, pure-CI and pure-orbital
directions through the orbital adjoint pullback.

Using only saved old measurements, the phase-instrumented LOFLEA run changes
165 sweeps to 129, 1639 generated batches to 1171 and 1,420,939,020 source
elements to 1,120,445,436. The complete coupled-image scope falls from
11.9900 to 9.7682 seconds; complete application time is 29.0093 seconds in the
native run. Four steps, 36 HVPs and the final energy are unchanged. The
standalone `hvp_apply` scope now includes the shared CI contraction; it must not
be compared with the old orbital-only scope as if their boundaries were equal.
See `results/diagnostics/c12-temporal-responses-20261017/` for raw evidence and
the saved old/native comparison. No old executable was rerun.
