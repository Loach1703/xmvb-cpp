# Orbital preconditioner and Armijo regression

Run only on Hanhai25 compute nodes. These benchmarks compare unchanged
mathematical preconditioners, not alternative optimization models.

`build.sbatch UPDATED_SOURCE` builds the executable and runs projector,
normalized-curvature, orbital-chart, and NEO globalization tests. The projector
test also has a `--benchmark` mode: 469 basis functions, 89 inactive orbitals,
11 active orbitals. It compares separate retained-span factorizations with one
shared factorization and full-AO rank-one corrections. This is a synthetic
preconditioner-only benchmark, not a molecular calculation or total speedup.

Submit `run.sbatch UPDATED_SOURCE BASELINE_SOURCE` with an `afterok` dependency
on that build. Its eight array tasks run in parallel: 241 (HAO) and FeCl2
(OEO), each with ISCF=8 and ISCF=9, using both binaries. Input files are the
same except for the requested optimizer. Both use projected-gradient tolerance
1e-3 and energy-change tolerance 1e-7 hartree, 32 OpenMP threads, and one BLAS
thread. Preserve baseline build flags and dependency versions.

Record source revisions/patches at submission; each task records input/binary
hashes, node, CPU, environment, raw output, exit status, and GNU-time wall time
and process peak RSS. Compare energy, convergence, outer iterations and action
counts before interpreting timings. A batch exit code alone is not convergence.

The implementation stores one full-AO projected model per chart, not one per
orbital. Sparse supports use the same factored inactive span without allocating
a full-AO model. Local curvature eigendecompositions remain unchanged.

## Phase and memory audit

Set `XMVB_PHASE_TIMING=1` to print aggregate CSV records to stderr on exit.
Diagnostics are off by default. `chart` includes projection, curvature and
block factorization; `block_factor` includes `block_eigen`. Do not sum these
inclusive totals. Projection totals include common setup and per-orbital models.
Trial totals include all objective trial calls, not just Armijo trials.
The projector workspace record measures only owned matrix payload, not RSS.

`phase-build.sbatch SOURCE_ROOT` builds `SOURCE_ROOT/baseline` and
`SOURCE_ROOT/updated` in parallel with the same Release flags, debug symbols,
and frame pointers. Both source trees must have the same timing instrumentation.
For the pre-optimization baseline, instrument the standalone projected-model
function rather than the new projector class; the measured operation is the same.

`phase-run.sbatch SOURCE_ROOT` runs four independent pairs in parallel. Each
pair evaluates baseline then updated on one exclusive node with the same CPU
allocation. Native time/RSS runs are separate from heaptrack captures; the
capture cases are 241/block-LBFGS and FeCl2/NEO. `phase-summary.py RUN_ROOT`
prints the paired resources and phases as CSV.

`phase-arena.sbatch SOURCE_ROOT` is an allocator ablation with
`MALLOC_ARENA_MAX=1`, not a proposed application default. It can substantially
slow multithreaded allocations and is not a production speed measurement.
`phase-memory.sbatch CAPTURE_ROOT` generates exact-byte Massif output and peak
allocation stacks from the saved heaptrack captures, on compute nodes.
`memory-summary.py CAPTURE_ROOT` extracts exact peak heap bytes. Chart bytes at
the global heap peak are not the chart's own maximum over its lifetime.
