# Shared sigma / integral-adjoint sweep

## Algebra

Let $X$ be the orthogonal spin-string product coefficient response and $Y$ the
accepted state, with RHSs packed in beta-column blocks. Let $D_t$ denote the
linear one-body scatter restricted to tile $t$. Its rows enumerate the tile's
spin-string products and RHSs; columns enumerate packed orbital pairs. The
absorbed integral kernel $K$ contains both one- and two-electron terms, using
the existing direct-CI normalization.

$$
\sigma(X)=\sum_t D_t^*\!\left[D_t(X)K^T\right],
\qquad
G_K(X,Y)=\sum_t D_t(X)^T D_t(Y).
$$

Here $D_t^*$ is the transpose scatter under the Frobenius inner product and
$G_K$ is the kernel adjoint, since

$$
\delta\langle X,H Y\rangle
=\langle G_K(X,Y),\delta K\rangle_F.
$$

The existing adjoint of the linear integral-absorption map converts $G_K$ into
one-electron and pair-kernel derivatives. Its symmetric Frobenius convention,
state weights, factor of two in response gradients, and energy-shift terms
are unchanged. The fusion changes evaluation order, not the physical operator.

Separate calls construct $D_t(X)$ twice and $D_t(Y)$ once. The fused sweep
constructs each once, accumulates the adjoint, then overwrites the consumed
$D_t(Y)$ storage with the sigma intermediate. Both required GEMMs remain.
This removes one scatter out of three; it does not lower the leading arithmetic
order of the integral contractions or imply a one-third wall-time reduction.

## Storage and ownership

The two shared scatter arrays together retain the existing 64 MiB payload
budget. Inputs, the required sigma output, the orbital-pair adjoint output,
and BLAS packing workspaces are separate from that payload; the budget is not
a bound on process RSS. No thread-private full coefficient carrier is added.

After the adjoint consumes its operands, the second tile changes layout in
place: each beta string owns a contiguous matrix. GEMM and alpha-link gather
consume that block consecutively. Beta links are gathered by output-column
ownership using the Hermitian packed one-body connection graph. Reverse links
have identical signs and pair indices. Barriers separate adjoint reads, reuse
of its storage, and the two output-writing phases.

Only a tile's live rows are consumed. This covers partial beta tiles and RHS
blocks without retaining a global scatter tensor. Scanning reverse beta links
does add per-tile index tests; this tradeoff must be measured.

## Interfaces and scope

`apply_with_adjoint(response, state)` returns sigma and integral adjoints.
`integral_adjoint` and this fused interface share `contract_adjoint`; the
standalone sigma remains necessary for Davidson actions that need no adjoint.
There is no retained duplicate adjoint implementation or runtime fallback.

`StructureAction::active_integral_response_adjoint` now uses the fused result.
Pure overlap multipliers still bypass response work by the existing exact-zero
rule. No optimizer settings, CI accuracy, HVP formula or convergence condition
was changed. Other consumers have not been claimed as fused.

## Verification protocol

All compilation and numerical work uses Hanhai22 Slurm, 32 bound OpenMP cores
and one BLAS thread. Direct-CI correctness tests also run at one thread.
Tests compare fused and separate operations for unequal vectors, self pairs,
zero response, multiple RHSs, non-contiguous input views, distinct alpha/beta
spaces, partial tiles, dynamic teams and nested OpenMP. Existing bilinear and
directional-derivative tests remain unchanged.

The bounded kernel benchmark compares both legitimate APIs in the same new
executable on identical synthetic carriers. It is not a rerun of an old build
or a full chemistry workload. Native convergence comparisons use saved data
from array 1250785, with identical input/basis hashes and thresholds.

The initial shared-column-layout prototype (build 1250803, kernel job 1250816)
passed correctness but its fused kernels were about 9–14% slower at 462/924
strings. The layout was replaced rather than kept as an alternative path.
Native array 1250806 showed no substantial end-to-end gain. Those artifacts
are retained as the unsuccessful experiment, not a performance baseline win.

The per-beta contiguous layout alone (build 1250829, kernel 1250830) did not
resolve the slowdown. Adjoint-only tiling still made a single beta task too
large and, for 924 strings with four RHSs, provided only 14 beta tasks for 32
workers. The retained implementation uses sigma's existing alpha extent of
160; this constant is shared with standalone sigma, not newly tuned. The
remaining shared budget holds more beta tasks. No rejected layout remains.

## Distinct spin sectors

Alpha and beta have independent counts and connection graphs. Links are shared
only when the complete occupation lists compare equal, never merely because
their sizes agree. Fixed spin sectors and the subsequent selected-VB-structure
projection are unchanged; no assumption of a full VB structure space is added.

Build 1250864 passed the eight test launches, including new cases at 1/32 threads:

| Active orbitals | Alpha/beta electrons | Alpha/beta strings |
|---:|---:|---:|
| 11 | 5 / 6 | 462 / 462 (different occupations) |
| 10 | 4 / 5 and 5 / 4 | 210 / 252 and 252 / 210 |
| 11 | 1 / 5 and 5 / 1 | 11 / 462 and 462 / 11 |
| 10 | 5 / 0 and 0 / 5 | 252 / 1 and 1 / 252 |

The extreme cases test the empty spin sector as well as transposed imbalance.
Small unequal-spin cases retain an independent explicit Hamiltonian reference.
This verifies correctness, not uniform parallel efficiency for all shapes.

## Kernel results

Final kernel job 1250869 compares separate and fused calls on the same new
binary, five calls per entry. Times below are milliseconds per combined call,
not end-to-end timings. The sector order in the archived CSV follows the
literal sector list in `benchmark_sigma`; two 462/462 entries differ in electron
counts even though their carrier dimensions match.

| Active orbitals; alpha/beta electrons | Strings | RHSs | Separate / ms | Fused / ms | Speedup |
|---|---:|---:|---:|---:|---:|
| 11; 5/5 | 462/462 | 1 | 36.797 | 35.035 | 1.05 |
| 11; 5/5 | 462/462 | 4 | 106.530 | 95.408 | 1.12 |
| 11; 5/6 | 462/462 | 1 | 26.783 | 27.925 | 0.96 |
| 11; 5/6 | 462/462 | 4 | 127.144 | 97.862 | 1.30 |
| 12; 6/6 | 924/924 | 1 | 119.958 | 134.691 | 0.89 |
| 12; 6/6 | 924/924 | 4 | 741.892 | 551.753 | 1.34 |
| 11; 1/5 | 11/462 | 1 | 0.808 | 0.717 | 1.13 |
| 11; 1/5 | 11/462 | 4 | 2.703 | 2.603 | 1.04 |
| 11; 5/1 | 462/11 | 1 | 2.191 | 1.378 | 1.59 |
| 11; 5/1 | 462/11 | 4 | 6.112 | 5.346 | 1.14 |

The fusion is not universally faster: the single-RHS 924/924 kernel is 12%
slower. Four-RHS improvements do not imply that the current per-state production
consumer benefits by the same amount. Do not hide these limits behind a
molecule-specific switch or report a blanket speedup.

## Native convergence results

Array 1250866, compared with saved 1250785 data. All six cases converged. Values
are before / after; no matched recent MnF2 baseline was available in that array.

| Case | Steps | HVP | Wall / s | Peak RSS / MiB | Final energy / Eh |
|---|---:|---:|---:|---:|---:|
| 241 | 3 / 3 | 21 / 21 | 2.069 / 2.049 | 820.7 / 797.9 | -230.720590392829 |
| FeCl2 | 12 / 12 | 256 / 254 | 7.963 / 6.700 | 491.2 / 503.7 | -2181.624464870615 |
| QOCGUM | 4 / 4 | 52 / 52 | 13.165 / 13.029 | 1915.4 / 2014.8 | -480.094766681857 |
| CERRAS | 5 / 5 | 51 / 51 | 18.174 / 19.974 | 1957.1 / 2012.3 | -397.093435004109 |
| LOFLEA | 4 / 4 | 34 / 34 | 34.131 / 33.999 | 2874.6 / 2821.5 | -399.106646001531 |
| MnF2 | — / 13 | — / 189 | — / 4.948 | — / 553.4 | -1348.893353223861 |

Maximum paired final-energy difference is 6.85e-10 Eh. Input/basis hashes match.
Node allocation differs for four paired cases, including FeCl2 and CERRAS, so
the FeCl2 timing decrease cannot be attributed wholly to fusion. CERRAS regressed
in this sample. There is no demonstrated across-the-board time or RSS improvement.
`response-comparison.csv` includes node names and unrounded values.

Diagnostic 1250870 (LOFLEA, two macro steps, 13 HVP expansions) confirms 12 fused
calls. Separate sigma calls decrease from 366 to 354; standalone adjoints fall
from 32 to 20. The fused scope takes 1.507 s, standalone adjoints 1.773 s, and
remaining sigma calls 9.057 s. These are diagnostic scopes, not production speed
ratios; sigma remains dominated by other consumers. The iteration-limit exit in
this diagnostic is expected; native convergence was verified independently.
