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
