# C12: one batch for retained directional panel reads

## Exact lifetime and index argument

An accepted retained state stores two numerical components on N supported
products plus its shared CSR metadata, and fits within the existing budget B.
Its directional input uses one component on the same support. Therefore
8N+metadata <= 16N+metadata <= B, even with conservative repeated metadata
charging. All directional input tiles fit without eviction. Gather only their
owner pointers before one collective beta-column copy; the pointers survive
every subsequent owner lookup in that read. Multi-root states obey this bound
individually. Other transient inputs retain panel-by-panel consumption, with
no pointer held across a possible eviction.

Sorted distinct alpha indices of length U_alpha are the identity row range.
Skip building/resetting its position map; full columns use contiguous Eigen
copies and compact columns scatter by their already exact original row index.
There is no amplitude threshold, new coefficient workspace, or change to the
support, physical HVP, CI response, optimizer, or convergence criteria.

## Qualification

Hanhai22 job **1264610**, `anode40`, 32 physical cores, affinity `60-91`.
GCC14.1, C++23 Release, original optimization flags and BLAS threads=1.
Completed exit0; **9/9 regression tests passed** in 3.43s. Existing narrow
workspace, selected-space, multistate and unequal-spin tests remain active.
Independent CAS12 errors:

| Check | Norm |
|---|---:|
| Selected H | 1.97942125e-14 |
| Selected M | 1.90486032e-15 |
| Physical structure-response image | 3.36632466e-16 |
| Structure-response integral adjoint | 1.47706187e-15 |
| Orbital HVP vs gradient finite difference | 5.84911563e-11 |

The optional complete-carrier product-jet replay is disabled; its placeholder
zeros are not used. The finite-difference checks execute independently.

## Converged native comparison

Original byte-identical LOFLEA CAS(12,12), COV132, HAO, exact cc-pVDZ, ISCF9,
Davidson and unchanged thresholds. Old versions are saved data, not rerun.

| Metric | Start: 1264598 | Gram: 1264607 | New: 1264610 |
|---|---:|---:|---:|
| Accepted steps | 4 | 4 | 4 |
| Orbital HVP / physical CI responses | 36 / 36 | 36 / 36 | 36 / 36 |
| Final energy / Eh | -399.106646267330 | -399.106646267331 | -399.106646267332 |
| SCF wall / s | 32.038107 | 29.927104 | 28.153246 |
| Application end-to-end / s | 34.899360 | 33.073310 | 31.147581 |
| Process peak RSS / MiB | 2697.812 | 2723.004 | 2687.047 |

The latest native run is 5.82% faster than the preceding saved run. Across
the three independently tested changes (accepted generator reuse, 2D Gram,
and panel reads), native wall drops **10.75%**, SCF wall **12.12%**. Peak RSS
is essentially unchanged versus the start; no new memory-scaling claim is made.
Every run takes 4 accepted steps and 6,17,10,3 orbital HVPs per step, terminating
with the original `ciah_projected_gradient_tolerance` condition.

## Separate phase-enabled results

| Inclusive scope / s | Gram: 1264607 | New: 1264610 |
|---|---:|---:|
| Application end-to-end | 33.071491 | 32.284880 |
| Complete coupled actions | 13.377332 | 12.851256 |
| Physical orbital HVP | 8.465292 | 8.011088 |
| Structure-side scope | 4.887558 | 4.819058 |
| Physical sigma/integral-adjoint | 10.910750 | 10.758141 |
| Tile reads | 3.695451 | 3.597307 |
| Tile compute | 6.364395 | 6.566289 |
| Gram, nested in compute | 2.954028 | 2.914878 |

The phase-enabled total improves 2.38%, less than the single native result.
Tile reads alone improve only 0.098s, so the full native difference must not
be attributed solely to reading. The strongest scoped gain across this whole
round remains 2D Gram: tile compute falls from 9.446894s at the starting
baseline to 6.566289s now. The saved direction-only phase run did not improve
uniformly and is retained in its own report rather than hidden.

All final phase runs retain 165 tiled sweeps / 3109 batches /
1,964,120,314 logical source elements and 114176 Gram output blocks. Nested
phase times are not disjoint percentages. No sampling-profiler or measured
DRAM-traffic claim is made. Remaining Gram and tile reads still cost about
2.91s and 3.60s, respectively; further optimization remains open.

`comparison.csv` includes all four new/saved native stages; `phase-summary.json`
and `evidence.tar.gz` retain parsed and raw logs, source/binary/input hashes,
environment, build and regression evidence. No unrelated RI changes were synced.
