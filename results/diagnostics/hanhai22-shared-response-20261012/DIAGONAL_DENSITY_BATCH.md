# Bounded exact-integral density batches for the analytic orbital diagonal

## Algebra and work reduction

The diagonal uses symmetric density sources, with

$$
\mathcal G[D]=2J[D]-K[D],\qquad
d_p(D)=\big[Q^T\mathcal G[D]Q\big]_{pp}.
$$

Here Q is the same accepted AO-metric orthonormal frame as before. Core
sources are the individual rank-one densities of its occupied core columns.
Active sources are the existing untruncated Householder-QR basis of the
one-/two-RDM density span. No singular-value/rank cutoff is introduced.
Linearity reconstructs every original active density channel from its QR
coordinates after contraction, exactly as in the previous scalar path.

All core and active sources now share interleaved exact-integral sweeps.
For m density sources and batch width b, integral decoding and shared-memory
graph reads decrease from m sweeps to ceil(m/b). Arithmetic remains
O(N_ERI m): this is batching/SIMD and storage optimization, not a new
asymptotic reduction of electron-integral contraction complexity.

## Packed output and memory bound

For symmetric physical AO response, workers accumulate raw opposite
orientations into a single packed AO-pair entry. Off-diagonal entries are
the sum of both orientations. Raw diagonal entries are doubled at final
reduction, preserving the original physical-output convention, including
coincident indices. Independent full four-index references test this rule.
Inputs, integral graph and result are shared; only bounded packed response
outputs and one channel scratch vector belong to individual workers.

Let N be the AO count, P=N(N+1)/2, t the worker count, and b the channel width.
The planned large temporary storage in bytes is

$$
B=8\{N^2+Nn_{\mathrm{active}}+
b[2N^2+t(P+1)]\}.
$$

It includes shared AO source/product scratch, the active contraction
temporary, interleaved source channels, full response columns, packed worker
outputs and their scalar-channel scratch. Small QR/source-coordinate arrays
and final MO diagonals are separate existing density-span data. Accepted
integrals are not replicated in this workspace or its workers.

The existing 64 MiB budget limits t and b. Width is also capped at the
CPU's double-precision SIMD lane count (AVX-512: 8, AVX2: 4); neither limit
depends on the molecule, gradient, energy or empirical convergence tuning.
If an exceptionally large AO basis makes even one channel/one worker exceed
64 MiB, that unavoidable quadratic working set is the lower bound and no
additional replicas are admitted. For the tested systems the full batch
workspace fits within 64 MiB.

FeCl2: 42 sources, width 8, 32 workers, 20.60 MiB planned workspace,
6 sweeps rather than 42 per diagonal build. Iron porphine has 112 sources;
the resource planner predicts width 2, to be checked against the actual
four-step job events. Do not claim its speedup before those jobs finish.

The core pair-map allocation is released immediately after its AO action,
before projection/exchange consumers. Core Coulomb projection output rows
are independent; their packed AO-pair vectors also use a bounded worker
budget. The shared AO-by-core product matrix is never copied per worker.

## Preserve the fast HVP workload

Preliminary 1249593 incorrectly used packed scatter for the single-direction
HVP too. It saved diagonal time, but expansion HVP time grew from 2.08 to
3.08 s and total time grew. That integration was removed.

The final implementation uses one templated contraction/reduction algebra
with two necessary storage roles: existing full-layout single-direction
HVPs and packed, pre-interleaved large diagonal batches. It is not a
failure-triggered fallback or numerical selector. The HVP uses its original
scalar stencil and summation order; the diagonal uses the packed channels.
RI keeps its existing signed low-rank density contractions and is not
converted to this exact-integral batching algorithm.

## Validation and native FeCl2 results

1249596: all five unit executables passed. Final build-only 1249601 passed
all five again, including the added two-core diagonal fixture.
Independent AO operator tests cover widths 1, 2, 3, 8, one/four workers,
repeated AO indices, packed/full output agreement, and packed worker shapes.
Existing diagonal tests compare every entry with independently transformed
four-index exact/RI references, including zero/dependent density spans,
signed weights and nonorthogonal accepted active frames.

Both native runs below are anode31, bound 32 cores, original guess, exact
integrals, identical HVP/CIAH/coordinate policies and convergence tolerances.
1249588 is archived evidence, not a rerun of the old implementation.

| Native run | Macros | HVPs | Energy / Eh | Wall / s | GNU peak RSS / MiB | Diagonal builds / s |
|---|---:|---:|---:|---:|---:|---:|
| Scalar sources, 1249588 | 12 | 249 | -2181.624464871953 | 8.43 | 503.77 | 1.4835 |
| Shared batches, 1249596 | 12 | 251 | -2181.624464871899 | 8.30 | 516.23 | 0.9065 |

Both have 13 diagonal builds. Measured build time falls 38.9%; total time
falls only 1.5%, so no large end-to-end gain is claimed. Peak RSS rises by
12.46 MiB while planned batch scratch remains bounded. Expansion HVP time
is 2.27 s after restoring its full layout, rather than 3.08 s in the rejected
packed-HVP integration. The final energy is 1.1834e-8 Eh above the PySCF
coupled-Newton reference -2181.624464883733 Eh. No accuracy condition changed.
This is a kernel-cost improvement, not a claim that the FeCl2 tail is cured.

Final numerical source is frozen in build 1249601. Regression array 1249602
(F2/241/MnF2) and bounded four-step iron array 1249603 use that binary and
have been submitted in parallel. Their results are pending scheduler/output
verification; they are not convergence claims. No old-version tasks are run.
