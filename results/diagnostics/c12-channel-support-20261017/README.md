# C12: exact channel-support audit, rejected local kernels

## Retained state

**No new production performance path is retained in this round.** The numerical
implementation is restored to d54c9095: `sigma_tiles.cpp`, `pair_index.cpp` and
the production CMake source list match that commit. The index header has a
documentation-only correction: its arguments/results are zero based, and pair
exchange is not antisymmetrization. Negative-argument checks are unchanged.

Retained additions are the bounded channel census, its small generator oracle,
the direct covalent-boundary channel certificate,
the derivation in `article/channel_support.md`, extra packed-index unit checks,
and inclusion of the existing common-integral regression in the C12 harness.
The three unrelated dirty local RI/ignore files were not transferred or staged.
Rejected implementations appear only inside the evidence archive, not in a
compiled production or maintained fallback path.

## New mathematical evidence

The ordinary generated-product support is not the support of every one-body
channel at a node. A channel can contribute only if its one-body preimage lies
in the certified source support. This is an exact graph certificate, not an
amplitude cutoff. Census jobs **1264972** and **1264993** streamed CAS(12,12)
even-exchange representatives, with **78** packed orbital-pair channels:

| Source family | Nodes | Mean live channels | Full-output multiply work fraction | Gram work fraction |
|---|---:|---:|---:|---:|
| One disjoint-pair covalent structure | 261868 | 18.5752 | 0.23814 | 0.11333 |
| Complete covalent family | 304590 | 24.8429 | 0.31850 | 0.17335 |
| FULL | 427350 | 61.3427 | 0.78644 | 0.61951 |

Fractions count mathematical multiply-adds, **not wall-time speedups**.
The single structure has 64 original source products; it is not the complete
covalent family. No full product mask is stored. Only one 4096- or 32768-row
diagnostic tile is retained, plus fixed-size global histograms and one-spin
maps. Final census 1265005 used **7764 KiB** process peak RSS and **16.29 s**;
these are diagnostic resources, not molecular SCF resources.

At 32768 rows per diagnostic tile, equal-mask groups average 3.37 / 1.99 / 1.17
rows. COV has 139318 singleton groups out of 153208 groups. Thus invoking a
tiny GEMM for each distinct mask is not justified by channel sparsity alone.
All grouping histograms, including the large groups, remain in the raw data.

Source/boundary separation is more informative:

| Family | Source nodes / mean channels | Boundary nodes / mean channels |
|---|---:|---:|
| Single | 68253 / 48.8445 | 193615 / 7.9047 |
| COV | 113322 / 50.9439 | 191268 / 9.3787 |
| FULL | 427350 / 61.3427 | 0 |

COV's boundary is **62.8%** of its nodes but averages only **9.38/78** channels.
This supports an interior/boundary algorithm rather than replacing every GEMM
with the same sparse loop. The existing `active_support` label already names
that distinction. No molecule-specific threshold is needed to define it.

The derivation distinguishes a universal full sigma image from an exact
selected bilinear form. Restricting both channel axes requires an explicit
left-consumer support contract. It cannot be inferred from a memory budget.
The current product-jet generator requests produce upper cotangents, whose
triangular pullback supports this route; transposed or symmetrized generator
requests must not silently inherit it. A complete consumer audit and actual
efficient implementation are still required. No current physical response
has been dropped on the strength of the census.

## An efficient boundary certificate, now verified

For the complete covalent source family, compute the first and last deficient
combined prefixes. An electron move repairs the node if and only if it moves
from an occupied site beyond the last deficit into an empty site at or before
the first deficit. Enumerate the two spin-specific Cartesian products with
bit operations instead of querying source membership on every orbital pair.
The derivation gives the exact channel count by inclusion-exclusion and an
O(n+k) emitted-index construction, with no per-pair overlap or cofactor work.

Job **1265004** verified the interval certificate. Final **1265005** verifies
the bit-set Cartesian-product generator against the general preimage mask on
**all 191268 COV boundary nodes**, independently with both diagnostic tile sizes.
There are no mismatches. Small n=4 generator-basis oracles also cover unequal
spin populations and the empty-spin sector. **11/11 regressions passed in
4.17 s** in final job 1265005. The maintained census source and target hashes
match this job.

This is a qualified **symbolic channel generator**, not a claimed production
sigma/gradient/HVP speedup. Its specialization requires the complete covalent
source family; arbitrary selected structures still require their actual
support/preimage graph. The next implementation must preserve the explicit
left-projection contract and all physical derivative consumers.

## Numerical support oracle

At n=4, including unequal spins and an empty spin sector, apply the production
one-body generator to every supported source basis vector and collect its exact
nonzero channel support. It equals the combinatorial certificate. This bounded
oracle has at most 36 carrier entries; it is not the n=12 storage algorithm.
All census versions report `small_generator_oracle=PASS`. The final bit-set
certificate additionally matches the unequal-spin boundary cases.

## A lower FLOP count did not make naive SpMM faster

Job **1264978** used 32 cores on anode42. Synthetic amplitudes populate the
certified channel masks in bounded 32768-row tiles. A current BLAS/Eigen GEMM
and a diagnostic exact-zero row-wise SAXPY implementation emit **all 78 output
channels**. Both are warmed before the recorded paired products. Three shared
numeric slabs occupy at most **58.5 MiB**; no worker copies a carrier.

| Mask family | GEMM / s | Exact-zero SpMM / s | Relative error |
|---|---:|---:|---:|
| Single | 0.006669427 | 0.008917591 | 2.2500e-17 |
| COV | 0.007457430 | 0.012212567 | 1.9746e-17 |
| FULL | 0.008297595 | 0.024371570 | 1.7558e-17 |

These are sums of tile-product times, **not SCF times or real accepted-point
kernel traces**. They reject the naive implementation, not the algebraic zero
certificate. The SAXPY loop repeatedly updates output rows instead of using
GEMM register/block reuse. No hardware cache/DRAM counters were captured, so
the table does not quantify memory bandwidth or a cache-miss cause. The
prototype is archived as `spmm_prototype.cpp` and removed from the maintained
census and production code.

## Complete-convergence tests of other rejected candidates

All runs below use Hanhai22 **anode42 / affinity96-127**, 32 physical cores,
GCC14.1/C++23 production Release flags, OpenMP close/cores and one BLAS thread.
The saved d54c9095 native control was **not rerun**. Actual molecular convergence
is LOFLEA HAO, exact Libcint/cc-pVDZ, 150 basis functions, CAS(12,12)/COV132,
ISCF=9/Davidson; it is not a separate carbon-C12 molecular geometry.
Original orbitals, thresholds, ITMAX=2000 and response policy are unchanged.

| Variant / job | Steps / orbital HVP / response | Energy / Eh | Native app / s | Peak RSS / MiB |
|---|---:|---:|---:|---:|
| Saved retained / 1264965 | 4 / 36 / 36 | -399.106646267332 | 17.764194 | 2720.47 |
| Alpha SIMD clone attempt / 1264980 | 4 / 36 / 36 | -399.106646267329 | 18.652149 | 2652.47 |
| Masked int-offset attempt / 1264985 | 4 / 36 / 36 | -399.106646267330 | 18.255819 | 2709.32 |
| Header constexpr index / 1264994 | 4 / 36 / 36 | -399.106646267330 | 18.131653 | 2653.78 |

None demonstrated an end-to-end gain. These are single measurements, not a
confidence interval or proof that every isolated instruction change is slower.
They are insufficient grounds to replace the retained production implementation.
All converge with per-step counts 6,17,10,3 and the original termination reason.
RSS is GNU-time process peak; it is not an attributed allocation measurement.

Separate phase runs are preserved in `phase-summary.json`. In particular, the
masked attempt reduced a diagnostic HVP scope but did not improve native total
time; that local result is not substituted for the end-to-end comparison.
The original two alpha attempts passed 9/9 tests; constexpr passed the enlarged
10/10 suite. Their source snapshots/diff are archived, not compiled defaults.

## Generated-code evidence

Report jobs **1264981** and **1264986** disassembled the frozen candidate binaries
on compute nodes. Neither intended gather loop contains packed-double FMA or
gather instructions; its floating arithmetic remains scalar. ISA dispatch and
an OpenMP SIMD pragma alone did not provide the intended vectorization.

Compute job **1264992** compiled a separate diagnostic object with the captured
production flags plus only `-fopt-info-vec-all`. Its report identifies unsupported
control flow and an indirect-double-load vectorization failure. A report line
saying another reduction/control loop used 64-byte vectors is **not** proof
that the hot gather vectorized. These observations do not isolate one sole
cause or prove that changing global spin-link layout would necessarily help.
The diagnostic object was never linked into a measured executable.

After removing all three production candidates, **tests-only job 1264997**
rebuilt the restored implementation: **10/10 passed in 3.56 s**. This job ran
no native SCF benchmark. The new packed-index ordering/error tests and the
existing exact/RI transform, adjoint and HVP tests passed together.
The subsequent final certificate job 1265005 expanded this to 11/11 tests;
it ran the new census, not an old native SCF timing control.

All installation/build/workload/report work was on Slurm compute nodes;
local work was source editing and saved-result parsing. Raw native/phase logs,
input/binary/source hashes, compiler flags and assembly are in `evidence.tar.gz`.
Executable and object binaries are excluded; frozen candidate executables
remain in their remote `c12-regression-JOB/bin/` directories. No final speedup
or memory reduction is claimed for this diagnostic-only round.
