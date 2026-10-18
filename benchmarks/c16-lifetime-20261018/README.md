# C16 column lifetime and shared workspace

This benchmark measures **only the first structure-space CI solve**, not a
converged orbital optimization. The external driver stops the process after
`ciah_ci_solved`; signal 15 and wrapper status143 are intentional. A timeout
is retained as a timeout, never reported as convergence. Active CI energies
are not molecular total energies. No historical executable was rerun.

## Changes

- Retire real input/output columns only after their last exact generator use.
- Share immutable exterior factor rows between expansion and projection;
  assemble rows contiguously and contract exact-zero CSC factors.
- Stream targets once with persistent support cursors. Admit their actual
  shared neighborhoods, hold at most one pending target, and build each tile's
  maps only once. No worst-case-degree rejection or target regeneration.
- Default the H/S action to4 GiB of shared numerical residency/scratch,
  configurable by `--ci-workspace-mib`. Product order splits the allowance
  between input residency, output residency, sigma and transformation scratch.
  Source order uses its own single bounded factor workspace, not four
  fictitious product partitions. Allocation is demand-driven, never per core.
- Leave thresholds, integrals, CIAH strategy, input guesses and physical
  structure/spin spaces unchanged. Other derivative consumers retain their
  explicit budgets. Immutable geometry data, integrals, Davidson vectors and
  allocator storage are outside this tile allowance; measure process RSS.

The derivation and payload formula are in
`article/exterior_column_schedule.md`. Planner tests cover selected/full
support, both exchange parities, pending-target boundaries, exact retirement,
and sigma equivalence across small budgets. H/S, adjoint, orbital derivatives,
CI transport and coupled-action regression tests exercise the actual kernels.

## Frozen workload and environment

Original `C16_vb.xmi`: C16H18 singlet, CAS(16,16), COV, HAO,
Libcint/cc-pVDZ;330 Cartesian AOs,1,430 structures and12,870 spin strings per
spin. CLI explicitly overrides input ISCF=8 with `--optimizer-backend ciah`.
Input SHA256:
`044a228603ede8cf877b72889f71e65e5d7d7a3d4e1ec950fd75f609ddd0d448`.

Hanhai25, anode017, AMD EPYC9755,32 physical cores, affinity224-255,
OMP close/cores, BLAS one thread. GCC14.3.0, C++23 Release,
`-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`, Eigen BLAS/LAPACKE enabled.
Each run retains its binary hash, source manifest, CMake cache, environment,
input, GNU-time RSS and raw scalar/profiling logs.

## Existing controls from this round

| Job | New-code control | First CI wall/s | CI iterations/actions | Active CI energy/Eh | Peak RSS/GiB |
|---|---|---:|---:|---:|---:|
|292938|Lifetime + shared factor rows, dense panel consumption|703.99|14/14|-20.806938222167037|6.046|
|292946|Above + exact CSC consumption, old tile planner|713.55|14/14|-20.806938222166981|5.780|
|292998,4096 MiB|One-pass actual-halo planner + shared4 GiB budget|187.18|14/14|-20.806938222167162|8.544|

These are new intermediate variants from this round, not new runs of an old
repository revision. CSC alone did **not** improve end-to-end time. Its
reduced scalar contraction count must not be presented as whole-CI speedup.

The final4096 MiB measurement completes the first CI solve3.81x faster than
the existing292946 control (73.8% less wall time), with identical14 block
actions and an active-energy difference of1.81e-13 Eh. Peak RSS increases by
2.76 GiB. This measures the combined scheduling and residency change; it does
not isolate memory capacity alone. The lower-budget runs use the exact same
new executable and are still in progress. No claim of complete CIAH/SCF
convergence is made by these first-CI captures.

The earlier short run292879 stored a95.60s instrumented first sigma sweep,
but ended at120s without completing CI; it is not a matched full-CI control.

## Sampling evidence and follow-up

Job292957 sampled the exact frozen binary of292946 with gprofng for75s on
anode003 (also AMD EPYC9755). Reports contain resolved functions and no timer
warning. Team barriers and mutex paths dominate the CPU samples; supported
left transforms, sigma scatter and exterior rows are the next visible work.
These are CPU sample shares, **not removable wall-time fractions**. This
motivates reducing repeated tile preparation and synchronization rather than
claiming an equivalent speedup from the barrier percentage.

`run.sbatch` qualifies and freezes a build on a Slurm compute node.
`budgets.sbatch` then runs the **same new executable** at4096,1024 and256 MiB,
sequentially on the same node, and captures one bounded4096 MiB phase trace.
Native timing is separate from instrumented timings. The lower budgets have
a600s diagnostic cap and may remain resource/time-limit observations.

Build-only292963 rejected an obsolete512-byte column-only test budget after
sigma was included in the same allowance.292994 exposed a source-order
partitioning error: its single factor workspace needs the unsplit allowance.
Both failed qualifications stopped before any molecular run. Tests now use
legal shared budgets; direct LRU and small-budget sigma tests still force
eviction and tile boundaries independently.292996 qualifies the corrected
source-order accounting; all6 tests pass in6.39s. The frozen binary SHA256 is
`57204727c0add6b806109a7d519af5b3e38d1958abb92c26afc9e61bcced386e`.

## One outer-step measurement

`step.sbatch` reuses this qualified frozen binary and the original input,
with 32 cores, 4096 MiB and `--max-iterations 1`. It records scalar CIAH events,
accepted-step diagnostics, phase timings and GNU-time RSS. Its 30-minute
external cap bounds the diagnostic; a timeout is not a completed step. The
first-CI-only driver is not used, so initial gradients, HVPs, candidate CI
solves and the first accepted orbital update are included. Phase timing is
enabled, so this is an instrumented one-step measurement. Full SCF convergence
is outside its scope.

Submitted as Hanhai25 job293003 on anode003 (AMD EPYC9755). The earlier
first-CI timing was on anode017; retain that node distinction when comparing.

## Product-jet completed-output scheduling

The first one-step capture293003 completes initial CI in184.03s, then spends
57.96s in the output callback of its very first gradient sigma tile. This is
not a measured full-gradient time. Code inspection identifies a residency
mismatch: H/S uses the configured4 GiB allowance, whereas product-jet
input/output panels and response sigma still use separate old64 MiB limits.
Output eviction invokes generator cotangents on partial images and rereads
neighboring transformed inputs. The control is not rerun.

The new derivative path reserves accepted data across roots, partitions the
remaining configured workspace, and retires completed output panels after
all temporal response/forcing writers. Input lifetimes are deliberately not
shortened: the generator adjoint still needs their neighbors. Fully resident
images are accumulated before exchange projection and reverse consumption.
Accepted two-lane payload counts shared structural indices once and must fit
the subsequent coupled-direction input allowance. The mathematical argument
is in `article/exterior_column_schedule.md`.

Build293117 caught a size-type deduction error before any workload ran.
Build293119 passes6/6 tests in7.38s and freezes binary
`d2fe908fa2e678b5582e1c970498a005abf14b0b93df876296132de028c83e8e`.
New tests check temporal-output retirement and compare images/integral
adjoints with independent complete-carrier oracles. Build293164 adds an
end-to-end32 KiB workspace regression on nonclosed selected spaces, both
exchange parities, unequal spins and multiple roots; numerical production
sources are unchanged from293119. The32 KiB test exposed an OpenMP retirement
race: a worker could test the mutable released-column list while another
started the next tile, entering different conditional barriers. The hung
qualification293164 was cancelled. Serial and phase-instrumented diagnostic
runs passed, showing why the ordinary large-workspace tests had missed it.
Retirement now executes inside the final writer's existing single/barrier;
there is no extra conditional barrier. Final qualification293166 passes6/6
including the32 KiB regression in10.96s. Final binary SHA256:
`7316f3a8b5414e907df7d597d9522a77b2eb183c137357676aad58ac0d66a1b8`.

New32-core one-step jobs293141 (native, anode014) and293142 (phase timing,
anode039) use the same frozen binary and original input. Job293143 takes a
420s gprofng capture of that new binary; samples are for attribution, not
native timing. No thresholds, CIAH policy, physical spaces or old executable
runs are changed.293141/293142/293143 use the pre-barrier-fix variant and are
retained as intermediate measurements. Final native job293167 uses the
qualified barrier fix. These jobs are in progress; do not infer one-step
completion or full SCF convergence from their submission or batch status.

In293142 the accepted-state sigma sweep completes with175 batches:
prepare4.0091s, read4.3625s, compute6.9511s, write2.0599s. The17.3827s sum
excludes the final panel-consumption work and is not a complete gradient
time. The completed accepted product jet reports66 input panels,66 output
consumptions, zero input/output evictions and695317744 retained bytes. It
has reached the CIAH macro model; the old293003 capture was still inside its
initial gradient. These different endpoints do not define a full-step
speedup factor.

Subsequent orbital-direction jets in293142 also report66 generated inputs,
66 output consumptions, zero evictions and the same retained accepted payload.
Native293141 reaches actual gradient keyframes (about20s) and orbital HVPs
(about30s each). These are intermediate-variant per-operation observations,
not complete outer-step timings or convergence results.

The420s gprofng capture293143 has resolved symbols, no collector errors or
timer warnings, and the expected `experiment was not closed` note after its
bounded timeout. It covers initial CI plus early CIAH work. Scatter, channel
GEMM, gather and team waits are the main sampled paths; generator-adjoint
replay is no longer a leading sampled function. CPU sample fractions are not
wall-time fractions and cannot be compared as a matched percentage reduction
against the earlier75s initial-CI-only capture.
