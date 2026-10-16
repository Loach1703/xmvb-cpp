# Bounded exterior input/output jets: end-to-end CAS12 qualification

## Necessary operation, unnecessary contraction order

The spin-string Hamiltonian/overlap action, integral adjoint and selected-space
restriction remain necessary. Their repeated beta-row reconstruction does not.
For a physical image Z the required projection is

    W^T vec(T_alpha^T Z T_beta).

The former product-jet consumer sampled every touched product against every
transformed beta panel. Across all panels this costs O(|Omega| B); broad
support has |Omega|=A B, giving O(A B^2). The resident group accumulator
already avoided repeated structure-incidence scans when it fitted. This round
does not incorrectly attribute that older cost to rescanning all structures.

Accumulate the alpha-pulled contributions first. With P=T_alpha^T deltaZ +
deltaT_alpha^T Z and Q=T_alpha^T Z, one directional beta adjoint returns
P T_beta + Q deltaT_beta, followed by one W^T contraction. The beta work is
O(n_active^2 A B), not O(A B^2). Separate H/M and diagnostic targets obey the
same linear identity. No selected-space closure is assumed.

The input producer had the dual O(|Omega| B) reconstruction: for each beta
panel it multiplied the raw touched amplitudes by the corresponding exterior
rows. The final version instead expands X=W c once and applies the ordered
alpha/beta exterior jet to X and deltaX=W q. Generated sigma readers borrow
its columns directly. Accepted Y/residual panels are reused; a pure CI
direction does not rebuild accepted Y. Orbital directions retain the primal
needed by the ordered derivative only until that derivative has been formed.

The input work is O(nnz(W)+n_active^2 A B). Although a pure orbital tangent
can also be expressed as E(deltaR R^-1)Y, the implementation uses the existing
ordered exterior jet, avoiding a new inverse-conditioned cancellation in
nearly dependent orbital charts. A generator-only prototype was not retained.

Equations and assumptions are in Sections 10.3--10.4 of
`article/source_column_hamiltonian.tex`; its PDF was compiled and inspected.
Sigma/adjoint channel sweeps, alpha/beta parity factors, finite Ritz residuals,
optimizer policy and all numerical thresholds are unchanged.

## Memory and selected spaces

This is an exact contraction-order choice within the existing four-quarter
structure workspace, not a memory cutoff that disables orthogonal CI.
Broad input uses at most two shared carrier arrays; broad projection uses
one/two per ordinary/directional target. Each family must fit its original
64 MiB quarter. Narrow support and larger carriers retain the bounded panel
schedule, still computing the same selected-space operator.

For A=B=924, one array is 6.514 MiB. A directional input needs at most
13.028 MiB, and the two production H/M projection targets need 26.055 MiB.
They replace the previous input/projection payloads rather than extending
the budgets. Roots are processed sequentially. No carrier or integral array
is duplicated per worker; no persistent cache was added.

## Correctness

Final build/test job 1263145 passed all 12 tests on Hanhai22. Added tests use
three arbitrary, nonclosed selected columns with broad touched support,
two states, and both equal and unequal alpha/beta dimensions. They compare
the resident order against a 128-byte generated-panel order and independent
complete-carrier H/M references. Orbital H/M directional differences pass
at absolute error below 2e-7; response, generator and connection-jet comparisons
pass at about 1e-11 tolerances. Independent accepted Y/(H-E)Y panels exercise
accepted-state reuse. Existing tests also cover singular occupied minors,
parity signs, zero weights, tiny nonzero directions, CI transport and keyframes.

Output-only build 1263104 also passed 12/12. It is retained as an ablation,
not confused with the final input/output implementation. No extra final
stationarity check or relaxed stopping rule was added.

## Native molecular one-step results

LOFLEA CAS(12,12), original HAO guess, cc-pVDZ, libcint exact integrals,
NEO/Davidson, analytic orbital Hessian diagonal, gradient 1e-3 and energy
1e-7 Eh. The saved remote is 946934570c096ba5072923effa7fb063736731f1.
It was never rerun.

All cases use 8 actual OpenMP cores within one LLC, BLAS 1, 16 CPUs allocated
and 48 GiB requested. FULL has 226512 selected structures, COV 132.

| Space/version/job | Node | HVP/response | Final CI actions | SCF (s) | Program (s) | Process wall (s) | RSS (GiB) |
|---|---|---:|---:|---:|---:|---:|---:|
| FULL saved remote / 1262623 | anode50 | 8/8 | 39 | 74.965957 | 80.973591 | 93.77 | 3.373 |
| FULL previous current / 1263023 | anode50 | 8/8 | 27 | 63.012854 | 68.972155 | 89.46 | 3.435 |
| FULL output-only / 1263111 | anode50 | 8/8 | 27 | 57.644907 | 63.529981 | 80.84 | 3.370 |
| FULL final input/output / 1263147 | anode04 | 8/8 | 27 | 54.074424 | 60.027518 | 75.75 | 3.350 |
| COV saved remote / 1262622 | anode50 | 6/6 | 99 | 53.464025 | 57.809377 | 58.06 | 2.612 |
| COV previous current / 1263025 | anode50 | 6/6 | 99 | 46.285648 | 50.527074 | 50.79 | 2.500 |
| COV output-only / 1263113 | anode50 | 6/6 | 99 | 44.957258 | 49.204730 | 49.45 | 2.510 |
| COV final input/output / 1263149 | anode04 | 6/6 | 99 | 39.691606 | 43.319410 | 43.61 | 2.518 |

The same-node output-only FULL result reduces SCF/process time by 8.5%/9.6%
relative to the previous current, and its process wall is 13.8% below the
saved remote. Thus there is an end-to-end improvement before changing nodes.

Anode50 subsequently filled with other users' jobs, giving a next-day start
estimate. Pending final jobs were moved to anode04, also dual AMD EPYC 9654
96-core processors with the same 8-NUMA topology. An attempted anode10 target
could not start because it was owned by another user; no calculation ran there.
Final native FULL uses CPUs 80--87, scalar FULL 160--167, COV 176--183;
these three jobs ran concurrently. Output-only native/scalar jobs used
176--183 sequentially on anode50; COV subsequently used the same group.

Final observed process times are 19.2%/24.9% below the saved FULL/COV remote
times, but those final comparisons cross nodes. Do not assign their entire
difference to this code change. In particular, sparse COV does not enter the
new broad-support schedule, so its time variation is not evidence for a new
COV kernel speedup. There are no repeated-run averages. Program timers omit
some final report work; GNU-time process wall includes the roughly 158 MiB
FULL output. Both are reported explicitly.

## Energies and iteration scope

| Space | Final-run initial energy (Eh) | After one final-run step (Eh) | Previous-current after step (Eh) | Saved-remote after step (Eh) |
|---|---:|---:|---:|---:|
| FULL | -398.880685140388 | -399.322366416767 | -399.322366416768 | -399.322429054264 |
| COV | -398.533053417516 | -399.074297760065 | -399.074297751296 | -399.074297757091 |

FULL retains 8 HVPs, 13 microsteps and 3 keyframes; COV 6, 10 and 3.
The FULL step norm is 0.3599442862314564, versus 0.3599442862314833
previously, and its one-step energy differs by about 1e-12 Eh.
The final physical CI residual is 6.34373739e-8 Eh, inside the unchanged
1e-7 Eh budget (previously about 6.34373535e-8 Eh).
The pre-existing +0.06264 mEh FULL offset from the remote one-step endpoint
is unchanged; see `full_ci_isometry.md` for the older coordinate-dependent
CI stopping change (39 versus 27 CI actions). This round does not take
credit for that earlier action reduction.

All jobs end deliberately at one accepted step with application exit 2 and
`max_iterations`; the Slurm application-step FAILED label is that limit,
not an SCF crash. These are end-to-end one-step benchmarks, not proofs of
converged-energy equivalence.

## Phase evidence and remaining work

| Eight actions, seconds | Previous / 1263024, anode50 | Output-only / 1263112, anode50 | Final / 1263148, anode04 |
|---|---:|---:|---:|
| Orbital HVP | 8.030054 | 4.654224 | 3.954305 |
| Integral-direction contraction | 7.190924 | 3.768061 | 2.953487 |
| Coupled structure phase | 4.247933 | 3.322648 | 3.459247 |
| Complete CIAH macro | 26.570133 | 22.304018 | 26.032298 |

The final scalar run has 63.396 s SCF time and 84.04 s process wall, versus
54.074/75.75 s natively. Its initial/final CI solves take 13.995/8.906 s;
the prior scalar run took 9.570/6.561 s with the same action counts. Thus
its slower complete macro must not be interpreted as a strategy regression
or hidden behind the faster HVP. Node/core-group contention and scalar
instrumentation are not isolated by these runs. Use native runs for runtime
claims, phase clocks to locate work, and do not add nested phase totals.
Earlier profiler-permission limitations remain; these are not sampled shares.

The beta reconstruction removed here is no longer the main coupled cost.
CI solves and the remaining sigma/adjoint channels are the next candidates;
their mathematical necessity and native costs should be checked before
another representation or cache change. The unchanged sparse COV schedule
needs its own attribution rather than extrapolating the FULL gain to it.

## Reproduction artifacts

`exterior_jet_scf.csv` is generated from saved application outputs;
`exterior_jet_phases.csv` contains machine-readable phase totals.
`exterior_jet_evidence.tar.gz` retains the scoped source patch, input/source/
binary hashes, build flags and 12-test results, native/scalar stderr, compact
stdout summaries and hashes, traces, CPU affinity and Slurm accounting.
Full native outputs remain locally and on Hanhai22; scalar FULL stdout
remains on the cluster. Unrelated local RI edits were not synchronized or
included in this exact-integral change.
