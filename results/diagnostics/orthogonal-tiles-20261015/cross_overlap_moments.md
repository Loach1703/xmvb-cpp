# Cross-overlap generator moments: necessary work and end-to-end qualification

## What the remaining 14 seconds actually contained

Current-only scalar timing job 1262975, built from a5ec68fd with timing events,
split the eight coupled structure actions as follows. The reference is a saved
current implementation, not a new run of GitHub 94693457.

| Eight actions, wall seconds | Before / 1262975 | After / 1263024 |
|---|---:|---:|
| Coupled H/M action and integral adjoint | 3.273610 | 3.229490 |
| Active-integral adjoint to orbitals | 0.019341 | 0.017925 |
| Coordinate projection and assembly | 0.017086 | 0.013567 |
| Residual-weighted CI transport curvature | 10.422144 | 0.984839 |
| Complete coupled structure phase | 13.734512 | 4.247933 |
| Complete CIAH macro phase | 41.868527 | 26.570133 |
| Two keyframes | 8.792220 | 4.198291 |

Transport curvature, not ordinary structure-coordinate projection or the
integral adjoint, was the dominant cost. The curvature contraction is 10.58x
faster; the complete coupled structure phase is 3.23x faster. Parent/child
timers overlap and must not be added. These are explicit application clocks,
not sampled-profiler percentages. The earlier VTune ptrace-policy failure
remains the applicable sampling limitation.

## Necessary mathematics; unnecessary reconstruction

For arbitrary selected structure expansion W, use the physical amplitudes
L=T_left W c_left and R=T_right W c_right. The cross-Gram action is exactly
W^T T_left^T [T_right W c]. The former source-column schedule reconstructed
all touched output products from transformed factors. For broad support this
included O(|Omega| r) reconstruction, with |Omega|=A B and r=B. It is
unnecessary: transform one shared amplitude and project once instead.

The orbital adjoint does not require that final structure projection at all.
Contract D_pq=<R|E_pq|L> directly using the existing alpha/beta one-body links.
Its chart adjoints are upper(D R_left^-T) and upper(D^T R_right^-T).
Differentiate the two amplitude arguments and the small inverse transpose to
obtain the exact directional adjoint. This includes the finite Ritz residual
in CI transport curvature; it is not set to zero or screened by magnitude.

The work becomes O(nnz(W)+n_active^2 A B), without a Hamiltonian sigma,
cofactor, explicit exterior matrix, or reconstruction through every source
column. This identity is valid for a nonclosed selected VB space. It does
not reinterpret the wave function as full CAS. Alpha and beta may have
different electron counts and string dimensions; combining their moments
requires only the same spatial orbital chart, not equal spin tables.

The formula and its directional chain rule are recorded in Section 10.2 of
`article/source_column_hamiltonian.tex`; the PDF was compiled and visually
inspected. The existing generator pullback shares the same small-matrix
chain rule through `moment_pullback()`.

## Memory and contraction order

There is no new persistent carrier cache and no per-thread carrier copy.
Broad touched support uses one/two shared carrier arrays for the primal/jet
action and two to four for the adjoint, only inside the existing caller's
64 MiB workspace budget. In this CAS12 case one carrier is 6.514 MiB and
the largest four-carrier payload is 26.055 MiB. Narrow support and larger
carriers retain the bounded source-panel contraction order. Independent
spin-specific charts also use that order. Both orders compute the same
selected-space operator; the budget changes scheduling, not the physics or
whether orthogonal CI is available.

Native FULL process peak RSS is 3.435 GiB, versus 3.354 GiB before this round
and 3.373 GiB in the saved GitHub run. Thus this change is a speed improvement
with a modest RSS increase, not a demonstrated reduction in peak memory.
The live carrier-workspace cap is unchanged; process RSS also includes other
live objects and allocator retention.

## Correctness qualification

Hanhai22 build/test job 1263019 passed all 12 selected tests. They include
independent occupied-minor cross-Gram references, directional-action and
adjoint finite differences, singular occupied minors, nonclosed arbitrary
three-column selections, unequal alpha/beta tables, CI transport, selected
metrics, coupled actions, shared Davidson and keyframe regressions.

The added five-orbital test compares the broad one-body-moment schedule with
bounded source panels: primal relative errors below 2e-12 and directional
adjoint/action errors below 3e-12. A 1280-byte primal budget and a 2560-byte
two-sided adjoint-jet budget force the independent source contraction order.
Zero and one-sided orbital directions are also checked. No precision,
optimizer decision rule, coordinate map, or numerical threshold changed.

## Native molecular one-step results

Original LOFLEA geometry/HAO guess, CAS(12,12), cc-pVDZ, libcint exact
integrals, NEO/Davidson, analytic orbital Hessian diagonal, gradient 1e-3
and energy 1e-7 Eh. All runs are on Hanhai22 anode50 with 8 actual OpenMP
cores within one LLC, BLAS 1, 16 CPUs allocated and 48 GiB requested.

| Space/version/job | Structures | HVP/response | Final CI actions | SCF (s) | Program timer (s) | Process wall (s) | RSS (GiB) |
|---|---:|---:|---:|---:|---:|---:|---:|
| FULL saved GitHub / 1262623 | 226512 | 8/8 | 39 | 74.965957 | 80.973591 | 93.77 | 3.373 |
| FULL prior current / 1262780 | 226512 | 8/8 | 27 | 78.151803 | 84.403653 | 104.20 | 3.354 |
| FULL new / 1263023 | 226512 | 8/8 | 27 | 63.012854 | 68.972155 | 89.46 | 3.435 |
| COV saved GitHub / 1262622 | 132 | 6/6 | 99 | 53.464025 | 57.809377 | 58.06 | 2.612 |
| COV prior current / 1262787 | 132 | 6/6 | 99 | 47.521203 | 51.664289 | 51.91 | 2.487 |
| COV new / 1263025 | 132 | 6/6 | 99 | 46.285648 | 50.527074 | 50.79 | 2.500 |

FULL improves over the prior current by 19.4% in SCF time and 14.1% in process
wall time, with identical action counts. It now beats the saved remote by
15.9% in SCF time and 4.6% in process wall. COV process wall is 12.5% below
the saved remote; its sparse support retains the source-panel order, so its
small timing difference from the prior current is not evidence for a kernel
speedup from the new broad-support schedule.

The program timer excludes some final output work; process wall includes it.
FULL stdout is about 158 MiB. Report both rather than using the shorter timer
as end-to-end process time. These are single shared-node measurements, not
repeat averages. New FULL native/scalar jobs ran concurrently on disjoint
CPU sets (80--87 and 176--183); COV subsequently used 176--183. The saved FULL
GitHub run used 176--183, while saved COV used 80--87.

### Energies and scope

| Space | Initial energy (Eh) | After one new step (Eh) | Prior-current after step (Eh) | Saved GitHub after step (Eh) |
|---|---:|---:|---:|---:|
| FULL | -398.880685140388 | -399.322366416768 | -399.322366416769 | -399.322429054264 |
| COV | -398.533053417517 | -399.074297751296 | -399.074297711215 | -399.074297757091 |

FULL retains 13 microsteps and 3 keyframes, COV 10 and 3. The new FULL
candidate norm is 0.3599442862314833, versus 0.3599442862292254 in the
before-timing run, and its final physical CI residual remains about 6.34e-8
Eh. The current-round end energy agrees with the prior current to 1e-12 Eh.
The larger 0.06264 mEh offset from remote is already present in a5ec68fd;
its cause and coordinate-dependent CI stopping behavior are documented in
`full_ci_isometry.md`. This round does not add a new tolerance change or
attribute the prior 39-to-27 CI action reduction to the moment contraction.

All molecular runs stop after one accepted macro step with application exit 2
and `max_iterations`. Slurm labels those application steps FAILED because
of exit 2; this is the deliberate one-step limit, not an SCF crash. These
tests qualify short trajectories and performance, not converged-energy parity.
The GitHub baseline was not rerun.

## Evidence

`cross_overlap_moments_scf.csv` is generated from saved outputs by
`tests/mcscf/summarize_scf_step.py`, selecting the six comparison rows above.
`cross_overlap_moments_phases.csv` contains the scalar breakdown.
`cross_overlap_moments_evidence.tar.gz` retains the scoped source patch,
source/binary/input hashes, build flags/log, 12-test output, native/phase
logs, compact FULL stdout extracts, input, trace, affinity and Slurm
accounting. Full native stdout is retained locally and on the cluster;
phase FULL stdout remains on the cluster. Unrelated local RI edits were
not part of this change or the exact-integral build.
