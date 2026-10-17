# C12: exclusive reusable scratch, no per-thread cache

## Necessary versus removable work

The coupled action's orbital and CI integral adjoints cannot simply be merged:
`ciah_step.cpp` passes their separate Hessian images to the existing coarse
preconditioner, for the independent directions (p,0) and (0,q). Their sum does
not determine those two images. Dropping one would change the optimizer and
is not a legitimate duplicate-work elimination under this task's constraints.

However, the channel scratch itself carries no accepted-point physical state.
Every input, channel, Gram and output view is written before use, and the API
already forbids retaining batch views across callbacks. A fresh allocation
and its OS first-touch are therefore not required for each sweep. The saved
parent's source allocates its planned slab in each nonempty `contract_tiled()`.

`ChannelWorkspace` now leases that same uninitialized storage. Its process-wide
pool retains at most **one idle slab**, capped at the existing **64 MiB** channel
workspace budget. No result is reused semantically, no per-thread numerical
cache or full carrier is introduced, and no physical work, precision, tile
geometry, convergence criterion or optimizer policy is changed.

The mutex covers only pointer checkout/return, never a numerical operation.
Concurrent or recursive calls obtain independent active buffers, as before.
Earlier returns free their buffer while any other lease is active; only the
last return may retain one bounded idle buffer. Thus an idle cached buffer
does not accumulate on top of a concurrent working set. An oversized buffer
is not retained. At times with no channel action, up to 64 MiB of scratch can
remain resident; this tradeoff is explicit and is not a claim of zero memory
cost. Ownership is swapped, not numerically copied. Fresh allocations remain
uninitialized using C++23's available `make_unique_for_overwrite` facility.

## Qualification

Hanhai22 Slurm job **1264887**, **anode42**, 32 physical cores, affinity
**96-127**, OpenMP close/cores, one BLAS thread, GCC14.1/C++23 Release production
flags. The old executable was **not rerun**: parent data are saved job 1264878
on the **same node and CPU set**. The unrelated local RI changes were not
transferred.

The operator fixtures are CAS(12,12)/COV132; full molecular convergence uses
**LOFLEA HAO**, exact Libcint/cc-pVDZ, 150 basis functions, 132 structures,
924 strings per spin, ISCF=9/Davidson. This is not a separate carbon-C12
geometry. Original input orbitals, thresholds and ITMAX=2000 remain unchanged.

Input SHA256:
`092d545498678b0cf772c5c84b9eab674b694375558cd391a1b739b7bdcc35eb`.
Basis SHA256:
`d53e7e02d55fb99cdf823d50cf13fec375e0cbe1220ff61816aee7ab38e38ec9`.

**9/9 regressions passed in 3.34 s**, including new simultaneous-call and
recursive-reader checks on the same immutable sigma action, different exchange
parities, multiple RHSs and repeated calls of different sizes. The dense
action/adjoint oracles and original tolerances are unchanged. Recursive and
concurrent tests check exclusivity and absence of deadlock, not parallel
scaling or a bound on allocator-internal retention.

A final tests-only job **1264894** strengthened the concurrency check with a
C++20 barrier in each call's first reader: both calls demonstrably hold live
scratch before proceeding, rather than merely being launched asynchronously.
**9/9 tests passed in 3.46 s**. Native source is unchanged; its hash matches
1264887 and the final test/production hashes match the tests-only manifest.
This job did not rerun a native SCF benchmark.

The CAS12 independent checks retain HVP/gradient finite-difference error
2.5933568576767072e-11, Lagrangian finite-difference error
6.9466100081384818e-12, structure-response image error
3.4002682088792616e-16 and integral-adjoint error 1.2707391900175511e-15.
The optional complete-carrier jet replay is disabled; its placeholder zeros
are not additional measurements.

## Native complete convergence

| Metric | Saved parent 1264878 | New 1264887 |
|---|---:|---:|
| Accepted steps | 4 | 4 |
| Orbital HVP / physical structure responses | 36 / 36 | 36 / 36 |
| Final energy / Eh | -399.106646267331 | -399.106646267331 |
| SCF wall / s | 16.871700 | 15.740349 |
| Application wall / s | 19.492537 | 18.280541 |
| GNU-time process peak RSS / MiB | 2719.13 | 2689.92 |
| Process user CPU / s | 505.24 | 477.33 |
| Process system CPU / s | 35.97 | 28.10 |
| Minor page faults | 4351525 | 2293315 |

Both use per-step HVP counts 6,17,10,3 and terminate with the original
`ciah_projected_gradient_tolerance`. The matched native application time
decreases **6.22%**, SCF time **6.71%**. Peak RSS does not increase in this
measurement; this is not a general asymptotic RSS-reduction claim. Single
matched measurements are not a statistical confidence interval.

Minor faults decrease about **47.3%**. They and process system CPU cover the
whole application, not just this allocator. Do not assign every fault or all
system time to scratch, or interpret minor faults as measured DRAM traffic.

## Separate diagnostic run

| Nested scope / s | Saved parent | New |
|---|---:|---:|
| Application | 19.333582 | 18.859224 |
| Coupled actions | 5.636303 | 5.397096 |
| Physical orbital HVP | 5.579662 | 5.339173 |
| Sigma/adjoint response | 3.612272 | 3.350102 |
| Ordinary sigma | 2.891516 | 2.566585 |
| Channel compute | 2.466888 | 2.464384 |
| Gram, inside compute | 0.410206 | 0.415089 |

Diagnostic counters show **1 allocation / 128 reused leases across 129
sweeps**, with maximum capacity **67108864 bytes**. Capacities are never summed
as if every reuse were another simultaneous allocation. The saved source path
would allocate once in each of these 129 nonempty sweeps; that old allocation
count is source/count inference, not an old profiler capture.

Counts are unchanged: **1171 batches**, **1120445436 logical source elements**,
**1068 product projections** and **25344 Gram output blocks**. Channel compute
is essentially unchanged, as expected: the change removes allocation/lifetime
overhead, not contractions or HVPs. The phase-enabled total gain is smaller
than the native gain and both are retained. Nested scopes overlap and must not
be summed as independent shares. No sampling-profiler percentage is claimed.

## Frozen anode40 historical follow-up

Earlier pending job **1264881** had not started and targeted the now-superseded
intermediate exchange-only candidate. It was cancelled **before execution**.
Replacement **1264892** completed on **anode40**, using `--bundle` with the exact
frozen **1264887** binaries, manifests and original input files. This is a
version update, not a timeout restart, and does not rerun an old baseline.
Binary hashes match across the two nodes; the delayed run did not pick up
later source changes. All CAS12 operator checks passed again.

| Metric | Saved 505e5405 / 1264871 | Frozen current / 1264892 |
|---|---:|---:|
| Node / CPU affinity | anode40 / 60-91 | anode40 / 0-31 |
| Steps / HVP / structure responses | 4 / 36 / 36 | 4 / 36 / 36 |
| Final energy / Eh | -399.106646267329 | -399.106646267333 |
| SCF / s | 19.845801 | 14.768960 |
| Application / s | 22.626016 | 17.208014 |
| Peak RSS / MiB | 2687.16 | 2633.84 |

Both are 32 physical cores on the same node, but the CPU/NUMA ranges and
background load differ. This table documents the historical end-to-end
observation, not an isolated percentage gain attributable to the algorithm.
Use the matched anode42 comparison above to evaluate the workspace change.
The anode40 phase run also records one allocation and 128 reused leases.

`comparison.csv` reuses the existing native parser. The common phase parser
now aggregates optional allocation/reuse counters and reports maximum capacity,
leaving missing older fields missing. Raw evidence contains source/binary/input
hashes, build cache, affinity, tests and native/phase logs, excluding binaries.

## Remaining open seam

Keep the two coarse-preconditioner images independent. A more promising next
dataflow target is the beta transpose gather: a fixed channel reads
`channel_values[target*pairs+pair]` with a packed-pair stride, inside the still
roughly 1.86 s output scope. Inspect whether an exact blocked output-layout
permutation can be produced during the existing GEMM copy-back, rather than
adding a third channel slab or a separate full transpose. This is an untested
next hypothesis, not a performance result or a license to change the model.
