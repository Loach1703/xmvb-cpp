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
