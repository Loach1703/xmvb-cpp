# C12 scatter addresses: retained and rejected candidates

Retained production change: [tile address plan](../../../article/tile_address_plan.md).
Relocate beta source addresses once per group, share them across alpha targets,
RHSs and response/state scatter passes, and omit identity full-alpha halo maps.
Keep the established single halo, row-major channels/images, GEMMs and barriers.
The physical arithmetic, optimizer, input, precision and memory budget do not
change. There is no full spin-product address table or thread-local replica.

## Verification

Hanhai22 job **1264636**, `anode40`, 32 physical cores, affinity `60-91`,
GCC14.1, Release C++23, original optimization flags and BLAS threads=1.
Completed exit0; **9/9 regression tests passed**, 3.35s. The compiled
production source and API match local hashes. Unrelated local RI edits were
not transferred. Independent selected H/M and physical-response references
and the CAS12 gradient/HVP finite-difference checks all pass unchanged.

Input hashes match the saved baseline: original LOFLEA CAS(12,12), COV132,
HAO, exact cc-pVDZ, ISCF9, Davidson. Runs complete convergence; old versions
are not rerun. Each version takes 4 steps, 36 HVPs (6,17,10,3) and 36 physical
CI responses and reaches the original projected-gradient termination.

| Native metric | Saved b1987c61 / 1264610 | Address plan / 1264636 |
|---|---:|---:|
| Final energy / Eh | -399.106646267332 | -399.106646267331 |
| SCF wall / s | 28.153246 | 28.207024 |
| Application wall / s | 31.147581 | 30.863387 |
| Process peak RSS / MiB | 2687.047 | 2824.676 |

The native wall decrease is only 0.91%; SCF is essentially flat (+0.19%).
Do not claim a large native speedup from this measurement. Separate phase
runs decrease application time from 32.284880s to 31.477191s; coupled actions
from 12.851256s to 12.414762s, orbital HVP from 8.011088s to
7.832145s, structure scope from 4.819058s to 4.561679s, and physical sigma/
adjoint from 10.758141s to 10.338633s. These inclusive scopes overlap.

Peak explicit tile payload is unchanged at **35.394608 MiB**. The compiled
beta map peaks at **11664 bytes**, shared by all workers. Complete alpha halos
cover **72.1004%** of measured generated targets. The extra integer map cannot
explain the native process peak increase of about 138MiB; no attribution to
allocator retention is asserted without a controlled test. The new phase-run
RSS is 2631.992MiB, demonstrating that process-peak variation must not be
confused with an explicit workspace reduction or increase.

## Rejected routes

| Candidate | Job | Native wall / s | Phase wall / s | Production status |
|---|---:|---:|---:|---|
| Shared response/state scatter, two already-budgeted halos | 1264619 | 32.428602 | 31.632713 | removed |
| Column-major Hamiltonian image | 1264632 | 31.751962 | 33.470527 | removed |
| Column-major image with contiguous weight runs | 1264635 | 32.989534 | 32.255964 | removed |

All passed the same 9 tests and retained 4 steps / 36 HVPs and the correct
energy. None beat the saved native baseline, so local timer improvements did
not justify retaining either production path. The last layout probe explicitly
measured scatter 1.270s, kernel/restore 2.820s, and weights 0.212s; its total
non-Gram compute and gather did not produce an end-to-end improvement.

No runtime switch, compatibility implementation, extra halo or new output
layout remains. Trial data/derivations and candidate diffs are archived only
as excluded experimental evidence. Final numerical sweeps remain 165, batches
3109 and logical source elements 1,964,120,314. No measured cache-traffic or
sampling-profiler claim is made from phase timers.

Remaining memory investigation: callbacks run in `omp single` sections, whose
owning thread can vary. A controlled primary-thread ownership experiment can
test allocator/workspace reuse without changing arithmetic. It must preserve
the existing barriers and task completion, and must be accepted or rejected
by native time/RSS rather than presumed faster.
