# Lifetime-exact halo planning: C12 convergence

Production change: `src/mcscf/ci/sigma_tiles.cpp`, explained in
[the derivation](../../../article/sigma_halo_budget.md). The planner previously
reserved two halos although their readers reuse one allocation sequentially.
Remove that obsolete liveness charge and derive the row ceiling from the real
two/three channel slabs. Keep all existing GEMMs, precision, optimizer policy,
64 MiB slab allocation and external panel-cache budgets. No new thread-local
carrier or persistent numerical cache is introduced.

Hanhai22 job **1264676** built and passed **9/9** regression tests (3.31s),
independent selected H/M and structure-response references, and CAS12 jet
finite differences. The HVP FD error is **5.8107905901e-11**; structure image
and integral adjoint relative errors are **3.4003e-16 / 1.6287e-15**.
Production source hash matches the compiled manifest. The added unit-test
code matches; its obsolete scratch comment was clarified after qualification.
Unrelated pre-existing local RI edits were not transferred or staged.

The matched molecule is **LOFLEA CAS(12,12)**, not a separately supplied carbon
C12 molecule: COV132, HAO, exact Libcint/cc-pVDZ, CIAH ISCF9 and Davidson,
original guess/ITMAX2000/thresholds. Input SHA256
`092d545498678b0cf772c5c84b9eab674b694375558cd391a1b739b7bdcc35eb`
and basis SHA256
`d53e7e02d55fb99cdf823d50cf13fec375e0cbe1220ff61816aee7ab38e38ec9`
match saved job1264636. Neither that nor any other old executable was rerun.

| Native metric | Saved1264636 | New1264676 | Same-new-binary1264685 |
|---|---:|---:|---:|
| Final energy / Eh | -399.106646267331 | -399.106646267332 | -399.106646267330 |
| Accepted steps | 4 | 4 | 4 |
| Physical HVP / structure response | 36 / 36 | 36 / 36 | 36 / 36 |
| SCF wall / s | 28.207024 | 28.794006 | 27.658785 |
| Application wall / s | 30.863387 | 31.677382 | 30.361555 |
| GNU-time process wall / s | 31.10 | 31.91 | 30.59 |
| GNU-time peak process RSS / MiB | 2824.676 | 2730.438 | 2775.211 |

All finish `ciah_projected_gradient_tolerance`, exit0; per-step HVP is
6,17,10,3. All use anode40, affinity60-91, 32 physical cores, GCC14.1 Release
C++23 and BLAS threads1. Repeat1264685 uses the frozen1264676 binary,
without rebuilding. This repeat resolves contradictory native/phase timings;
it does not justify selecting only the fastest new time. Native spread spans
the old result, so **no stable end-to-end speedup is established**.

## Separate phase-enabled evidence

| Quantity | Saved1264636 | New1264676 |
|---|---:|---:|
| Numerical sweeps | 165 | 165 |
| Batches | 3109 | 1639 |
| Logical source elements | 1964120314 | 1420939020 |
| Prepare / s | 1.558851 | 1.331151 |
| Halo read / s | 3.412612 | 2.644525 |
| Compute / s | 6.165370 | 6.374196 |
| Gather/write / s | 2.328785 | 2.129010 |
| Coupled images / s | 12.414762 | 11.990040 |
| Physical orbital HVP / s | 7.832145 | 7.569020 |
| Structure scope / s | 4.561679 | 4.396632 |
| Application / s | 31.477191 | 30.913943 |
| Process RSS / MiB | 2631.992 | 2795.367 |
| Explicit tile payload peak / MiB | 35.394608 | 60.944092 |

Timers overlap and are not sampling percentages. The unchanged allocated
slab remains64MiB; larger occupancy is not a larger budget. Fewer halo reads
and address plans are real count reductions, but GEMM work is unchanged and
compute time did not improve. The simpler planning proof is retained without
claiming a large total-time gain or a reliable process-RSS reduction.

## Removed experiment

The two-slab in-place prototype (jobs1264662/1264664) passed ten tests but
took31.462584/31.070410s native versus30.863387s saved. Its dead-halo scratch
helper, test target and production path were deleted, not retained as a
fallback. Its source diff and numerical evidence are in the sibling
`c12-two-slab-20261017` archive. A setup-only first job1264660 failed because
Make had not refreshed a new test target; explicit CMake configuration now
precedes the build. Job1264684 was canceled while pending to correct its
one-off report script; it ran no workload.

Raw logs, hashes, CMakeCache, tests, environment and inputs are archived without
executables. The frozen binaries remain remotely under
`/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015/c12-regression-1264676/bin`.
`comparison.csv` contains all native/phase candidate rows, not just the fastest.
No profiler-percentage or large-scale extrapolation claim is made.
