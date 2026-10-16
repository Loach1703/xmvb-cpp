# 240 FULL: current 32-core run against the saved excellent baseline

## Result

The current implementation is not slower than the saved equal-energy32-core
baseline. Current SCF time is3.082727 s versus3.397351 s (9.26 percent lower),
with the same4 accepted steps and30 orbital HVPs. Both print the final energy
-343.446302545202 Eh. This is the equal-energy historical reference, not a
claim to beat every earlier optimizer record.

No C++ algorithm, tolerance or optimizer option changed for this audit.
Current source revision is5618682f. Reuse of the previously qualified current
binary was checked against job1262005's binary and source manifests; its SHA256
is7b7a10d07f9402685c92d373ab10063b2249dcf21aba03cb8007742db3f0eb24.
No historical executable was rerun, and no local/login-node workload ran.

## Matched physical and execution settings

- System:240 FULL, HAO, CAS(8,8),1764 VB structures, one state.
- Exact Libcint integrals, cc-pVDZ; NEO/CIAH, structure-space Davidson and
  analytic Hessian-diagonal orbital preconditioning.
- Unchanged projected-gradient tolerance1e-3 and energy-change tolerance1e-7 Eh.
- Both input SHA256:60ba15a98b963caa1139676e661b2caacd0df9988a00d3fc74b71c8e268197d8.
- Both basis SHA256:d53e7e02d55fb99cdf823d50cf13fec375e0cbe1220ff61816aee7ab38e38ec9.
- Hanhai22,32 OpenMP workers, cores/close binding, maximum active levels1,
  BLAS threads1. Both hosts have AMD EPYC9654 96-core processors, two sockets,
  one hardware thread per core. Old hostanode62, new hostanode21. New allowed
  CPUs96--127 contain32 physical cores; the old allowed set is136--167.
- Current Release flags remain-O3 -DNDEBUG -g1 -fno-omit-frame-pointer;
  no profiling rebuild or new ISA flags. Shared-node single native measurements,
  not repeated statistical timings or an isolated-node strong-scaling study.

anode62 lacked32 free cores at submission. The initially requested anode14
had idle cores but an exclusive Owner belonging to another user. The same
pending job1262045 was moved to unowned anode21 before execution. No task ran
on anode14, no duplicate production run was submitted, and this scheduling
correction did not change requested cores or numerical settings.

## Native comparison

| Quantity | Saved1259844 | Current1262045 |
|---|---:|---:|
| OpenMP cores | 32 | 32 |
| Accepted steps | 4 | 4 |
| Orbital HVP | 30 | 30 |
| Initial energy / Eh | -343.075645280883 | -343.075645280883 |
| Final energy / Eh | -343.446302545202 | -343.446302545202 |
| SCF wall / s | 3.397351 | 3.082727 |
| Application end-to-end / s | 4.968763 | 4.536277 |
| Process wall, GNU time / s | 5.11 | 4.67 |
| Process peak RSS / MiB | 1641.074 | 1570.262 |
| User+system CPU / s | 119.21 | 106.98 |
| Mean accounted CPU cores | 23.33 | 22.91 |

Mean accounted cores include runtime spin/wait CPU and are not a measurement
of useful arithmetic occupancy. Slurm step MaxRSS is0 for this short job;
use native GNU-time high-water RSS, not that unpopulated accounting value.
The native peak is70.81 MiB lower, but no allocation-profile attribution to
a specific source change is inferred.

The earlier9.936526 s result used8 workers. Comparing it with this32-core
historical reference did not establish an implementation regression. The
immediate-parent8-core control remains useful for its own matched comparison,
but is not the primary historical performance target.

## Phase attribution

Separate current scalar-phase run:4 steps,30 HVP, same printed final energy,
3.104386 s SCF. Historical scalar-phase SCF is3.413957 s. Instrumented times
below are not substituted for native times. Nested rows are not additive.

| Phase | Old calls | New calls | Old wall / s | New wall / s |
|---|---:|---:|---:|---:|
| Orbital HVP | 30 | 30 | 0.523008 | 0.410006 |
| AO response, inside HVP | 30 | 30 | 0.356821 | 0.304118 |
| Selected direction | 30 | 30 | 0.049306 | 0.020913 |
| Active backward | 30 | 30 | 0.021642 | 0.009432 |
| CIAH coupled images | 30 | 30 | 0.657639 | 0.621251 |
| CIAH structure portion | 30 | 30 | 0.128463 | 0.201273 |
| CIAH keyframes | 10 | 10 | 0.676560 | 0.625007 |
| Trial forward | 5 | 5 | 0.608748 | 0.520642 |
| CI sigma | 228 | 183 | 0.131679 | 0.076959 |
| CI integral adjoint | 45 | 30 | 0.013597 | 0.003730 |
| CI generator adjoint | 518 | 150 | 0.033008 | 0.003656 |
| Transport metric solves | 32 | 32 | 0.022565 | 0.020479 |
| Physical metric actions | 190 | 190 | 0.015685 | 0.038194 |
| Physical metric solves | 9 | 9 | 0.017486 | 0.032792 |

The orbital HVP total is21.6 percent lower in these captures. Its AO response
and backward portions improve, as do the sigma/integral/generator contractions.
The per-outer-step orbital HVP counts remain7,5,15,3; keyframe/microstep counts
and rejected-trial pattern are also unchanged. This case therefore does not
show optimizer-strategy/action-count deterioration.

The largest local regression is the30-call CIAH structure portion:+0.072810 s,
or about2.43 ms extra per call. Its source boundary in
`src/vbscf/optimization/backends/ciah.cpp` covers the selected CI Hessian/metric
action, response integral adjoint, reverse orbital coupling and subsequent
coordinate handling. It is not the standalone sigma timer. The current
`StructureAction::contract_response` uses the generated panel sigma/adjoint
sweep rather than the old dense response expansion. Sparse incidence expansion,
column-panel ownership, generated contractions and projection are actual
operations in this chain. They are candidates for focused reuse/fusion work;
their individual contribution to the73 ms excess is not yet independently
timed or allocation-profiled. Do not claim that a particular map, allocation
or barrier accounts for the whole excess.

The `ci_response` counter increases30 ->45, but it does not mean15 new physical
outer responses. Old accepted-state preparation separately called sigma and
integral_adjoint; the new preparation uses fused apply_with_adjoint, sharing
the CiResponse timer with tiled coupled responses. The15-call counter transfer
is accompanied by integral_adjoint45 ->30. Actual coupled structure directions
remain30. Do not diagnose duplicate physics from changed timer classification.

Transport remains the factorized complete-selected-space inverse:32 calls,
about20 ms total, with no iterative selected metric solve. Thus removing the
optional restricted-completion preconditioner in5618682f is not the direct
cause of this FULL result. The two metric-norm counters also keep their call
counts; their wall increases are local costs, not new repeated certifications.

## Profiling limitation

All collection/probe/report commands ran inside the same Slurm compute step,
after native and phase runs. `perf` is absent from the loaded compute-node
PATH. Intel VTune2022.2 software hotspots was attempted, but refused before
launching the application because node ptrace scope is restricted. Capture
exit is1, stdout empty, and no VTune result exists. Report errors and collector
stderr are retained; batch completion is not reported as profiling success.
No system setting or privilege change was attempted. No sampled CPU hotspot
percentages are available. Conclusions here use actual phase boundaries and
counters, not a fabricated sampling profile.

## Remaining target and next decision

1. Keep both historical-excellent and immediate-parent comparisons mandatory;
   the gate is recorded in root TODO_ORTHOGONAL_CI.md.
2. Equal-energy FULL32 is already below the3.40 s reference. Do not describe
   the8-core9.94 s measurement as a32-core regression.
3. If optimizing this small FULL case further, inspect the measured73 ms
   structure-response excess. Preserve the bounded arbitrary-selected-space
   contraction; do not recover a small-case constant by reintroducing full
   product ownership or changing precision/optimizer thresholds.
4. Earlier Hanhai25 records near1.48 s SCF/12 HVP remain a separate historical
   target. Their endpoint differs by about2e-7 Eh and their optimizer trajectory
   differs. This audit does not establish superiority to that faster target.
5. This FULL control does not qualify large arbitrary selected spaces or
   complete the broad-support gradient/HVP memory migration.

`compare_full32.py` generates `full32_comparison.csv` and
`full32_comparison_phases.csv` directly from saved records. Current raw output,
phase logs, hashes, flags, affinity and failed-profiler evidence are archived
with Slurm accounting; the historical artifacts are reused without rerun.
