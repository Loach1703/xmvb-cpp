# C12 contraction audit: exact identity elimination and output batching

## Scope and reference

Hanhai22 Slurm, **anode42**, 32 physical cores, affinity **96-127** for every
native entry in `comparison.csv`. GCC14.1, C++23 Release `-O3 -DNDEBUG`,
debug symbols/frame pointers retained, OpenMP close/cores and single-threaded
BLAS. No local workload execution or old native baseline rerun was performed.
These are single measurements of different candidates, not repeated samples
or a statistical confidence interval.

The operator fixtures use CAS(12,12)/COV132. Complete molecular convergence is
**LOFLEA HAO**, exact Libcint/cc-pVDZ, 150 basis functions, 132 structures,
924 strings per spin, ISCF=9/Davidson. Here C12 names the existing 12-active-orbital
workflow, not a separate carbon-C12 molecular geometry. Original orbitals,
tolerances, ITMAX=2000 and physical response policy are unchanged.

Input SHA256:
`092d545498678b0cf772c5c84b9eab674b694375558cd391a1b739b7bdcc35eb`.
Basis SHA256:
`d53e7e02d55fb99cdf823d50cf13fec375e0cbe1220ff61816aee7ab38e38ec9`.

The controls are saved job **1264887**, revision **338389f2**, and saved parent
**1264878**, revision **d1d85b3e**. They were not rerun. Successful profiler
capture 1264928 sampled the then-current retained 338389f2 executable and is
not a new native timing control.

## What remains mathematically necessary

Orbital and CI integral adjoints cannot simply be summed and one discarded:
the unchanged coarse preconditioner consumes their separate images for the
independent directions (p,0) and (0,q). Their sum does not determine both
images. Common accepted-state channels and generated tile graphs are already
shared by temporal responses. The selected-structure space is unchanged;
no CAS closure is assumed, and no outer response or nonzero coefficient is
screened out.

## Retained implementation

1. Strictly zero primal shear coefficients represent exact identity factors
   and skip their primal forward/adjoint/inverse operations. **No small-value
   tolerance is introduced.** Directional transforms keep every original edge:
   a zero primal coefficient may have a nonzero derivative.
2. Product-jet output writers borrow a group of compact panels whose exact
   numerical/CSR payload fits the **existing** output-cache budget. They then
   launch the existing disjoint beta/RHS/alpha-stripe tasks once per group,
   rather than once per panel. Only a reusable small list of panel pointers
   and local column indices is new. No worker copies a panel. LRU entries
   touched by the batch cannot be evicted while its complete payload fits;
   tasks finish before the next group changes the cache. The original
   single-oversized-panel rule is unchanged.

Supported exterior transforms retain their existing column-local traversal;
shared-RHS alternatives are not maintained. No additional edge graph,
numerical slab or per-thread coefficient copy is retained.

The corresponding algebra and lifetime argument are in
`article/exterior_output_batches.md`. This is exact repeated-work and task-group
elimination, not a lower asymptotic physical-HVP scaling claim.

## Native complete convergence

| Metric | Saved 1264887 | Retained new 1264965 |
|---|---:|---:|
| Accepted steps | 4 | 4 |
| Orbital HVP / structure response | 36 / 36 | 36 / 36 |
| Final energy / Eh | -399.106646267331 | -399.106646267332 |
| SCF wall / s | 15.740349 | 15.287966 |
| Application wall / s | 18.280541 | 17.764194 |
| GNU-time process wall / s | 18.48 | 17.96 |
| GNU-time peak RSS / MiB | 2689.92 | 2720.47 |
| User / system CPU / s | 477.33 / 28.10 | 461.39 / 31.55 |

Observed application time decreases **2.82%**, SCF time **2.87%**. The process
peak RSS is **30.55 MiB higher**, not lower. The source introduces no large
numerical allocation, but this does not establish invariant RSS or explain
the measured peak difference; no allocation attribution was captured for this
candidate. Do not relabel that increase as a measured memory reduction.
Relative to the older saved pre-workspace parent, application time is
19.492537 -> 17.764194 s; this includes the earlier scratch-lease improvement.

All variants in the CSV converge in 4 steps with per-step HVP counts
6,17,10,3 and the original `ciah_projected_gradient_tolerance` termination.
The final-energy spread is at most 3e-12 Eh at printed precision. Optimizer
settings and thresholds were not changed to obtain these times.

## Separate phase-enabled diagnostic

| Nested scope / s | Saved 1264887 | Retained new 1264965 |
|---|---:|---:|
| Application | 18.859224 | 18.190408 |
| Physical orbital HVP | 5.339173 | 5.119607 |
| Sigma/adjoint response | 3.350102 | 3.120510 |
| Ordinary sigma | 2.566585 | 2.504473 |
| Channel output callbacks | 1.861408 | 1.750314 |
| Supported/full left transforms | 0.681449 | 0.670686 |

The physical-HVP and response scopes decrease about **4.1%** and **6.9%**
in these diagnostic measurements. Do not attribute all changes to a single
loop. Nested scopes overlap, and timings of unchanged
subkernels fluctuate too. The native result is the end-to-end comparison.

Counts remain **129 sweeps**, **1171 batches**, **1120445436 source elements**,
**1068 product projections**, and **1 scratch allocation / 128 reused leases**.
Maximum pooled channel capacity remains **67108864 bytes**. No full product
buffer is newly introduced. Shared-RHS-only job 1264950 recorded 18.142712 s
native. Adding output batching in 1264964 gave 17.690658 s; its left-transform
aggregate still took 0.740276 s. The final simpler 1264965 removes RHS loop
fusion: native time differs by only 0.074 s (less than 0.5%), while left-transform
diagnostic time is 0.670686 s and native peak RSS is lower. There is no clear
general advantage to retaining that extra loop transformation, so it is removed.
All these individual measurements, including the slightly faster intermediate
total, remain in the CSV; no minimum is silently selected as the final result.

## Numerical qualification

Final job **1264965**: **9/9 regressions passed in 3.34 s**. Existing tests
exercise unequal spin counts, exchange parities, arbitrary selected structures,
multiple states including zero weights, nonresident accepted data and **128-byte
output budgets forcing repeated eviction**. New checks compare a multi-RHS
supported transform with independent one-RHS actions and retain a nonzero
directional derivative at the identity transform using finite differences.
Temporal callback ordering, an adjoint-only middle request and a subsequent
reader exception are also covered.

CAS12 independent checks retain HVP finite-difference error
2.5933568576767072e-11, Lagrangian finite-difference error
6.9466100081384818e-12, response image error 3.4002682088792616e-16 and integral
adjoint error 1.2707391900175511e-15. The optional full-carrier jet replay is
disabled: its printed placeholder zeros are **not** extra measurements.

## Rejected candidates: no production alternatives retained

| Candidate / job | Native application / s | Decision |
|---|---:|---|
| GEMM orbit weighting / 1264911 | 18.686125 | Removed |
| Block-column output layout / 1264919 | 19.693073 | Removed |
| Fuse response writer and next reader / 1264943 | 18.698483 | Removed |
| Compact beta target iteration / 1264949 | 18.667509 | Removed |
| Per-shear integer edge list / 1264963 | 18.566209 | Removed |

All passed numerical regressions, but none demonstrated a native end-to-end
gain over the saved retained control. The blocked layout also increased the
output diagnostic from 1.8614 to 3.1690 s. Compact beta iteration is not an
automatic win: COV writers already use compact selected rows, so a second
ordered target list and intersection can add work instead of removing it.
The edge-list alternative adds metadata work without a demonstrated total gain.
Snapshots/derivations are scientific evidence only, not compiled or maintained
fallback implementations. `sigma_tiles.cpp` is exactly its retained parent;
the workspace lease remains intact.

## Sampling evidence and reproducibility

Initial capture-only job 1264920 could not use perf/VTune; the script's zero
exit status is not a successful profile. On compute node anode42, job 1264927
installed user-local GNU gprofng **2.42**, without sudo or kernel changes.
The source archive's SHA256 and build/install logs are preserved; no claim of
independent release-signature verification is made. Job **1264928** then
successfully sampled the frozen retained executable. Its version, binary/input
hashes and quality reports are archived. Report-only job 1264945 did not run
the application again.

The successful report records **no errors and no warnings**, a 10007-us
interval, 19.963629-s duration, and samples from **all 32 OpenMP threads**.
Sampled CPU sums to about 545.75 s versus resource-accounted 546.16 s. These
checks support using resolved hotspots, not claiming perfect unwinding or
turning CPU-sample shares into wall-time fractions. The 19.780960-s profiled
application is never used as a native control.

Resolved self samples include alpha transpose gather at sigma_tiles.cpp:682
(13.059 CPU s), beta gather at :703 (7.725 CPU s), and substantial AO row,
canonical-stream and AO-pair contractions. Alpha gather, not only the apparent
beta stride, therefore needs attention. OpenMP barriers have many CPU samples,
but that is neither an equivalent removable wall cost nor proof of imbalance.

Core `.er` sample data and all reports are in the evidence archive. The full
175-MiB capture, including archived binaries, remains at:

```
/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015/
  c12-cost-1264928/loflea-profile/capture.er/
```

`capture.sha256` inventories the full remote capture. Local copies omit its
`archives/` binary payload only; exact executable hashes and the frozen native
bundle remain available remotely. Candidate binaries are in each remote
`c12-regression-JOB/bin/`; they are not committed to Git.

The updated `tests/mcscf/c12_profile.sbatch` in `--capture-only` mode skips native operator
and SCF reruns, preserves the qualified bundle's build/source manifests, prefers
usable perf then user-local gprofng, and records function, line and thread
reports. Installation/probing, collection and report generation ran under Slurm;
only source work and saved-result parsing were local.
