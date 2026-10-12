# Matched CIAH preconditioner suite

## Protocol

Hanhai22 Slurm array **1247682** covers 15 inputs (13 HAO, 2 OEO), including
10 additional inputs beyond the previous five-system check. Cases are listed
in `preconditioner-suite.tsv`; all are single-state calculations. State-averaged
inputs are intentionally excluded because the current multistate NEO backend
is different from the single-state CIAH implementation being examined.

For each case, analytic Hessian-diagonal and positive one-electron are run
sequentially in one allocation, on the same node and 32 bound physical cores.
Cases are independent parallel array tasks; launch order is scheduler-dependent.
The order of the two preconditioners alternates across cases. Both use the same
frozen executable, exact integrals, Davidson eigensolver, orbital input and
NEO backend. No chemistry is executed locally or on a login node.

The one-electron variant is the restored **unshifted positive block inverse**,
not the problematic signed one-electron AH inverse and not the experimental
complete-metric shifted model. No optimizer code changes were made this round.

Frozen binary: `ciah-positive-build-1247629/xmvb-cpp.exe`, implementation
`ab7e3c09`, SHA256
`10430688c99848326f6c21967d4bf7e43c4a2ba1db1628b6e7a63b1e064e98a1`.
Compiler/build flags are retained alongside each paired run.

Threads: OpenMP 32, `OMP_PLACES=cores`, `OMP_PROC_BIND=close`,
`OMP_MAX_ACTIVE_LEVELS=1`; OpenBLAS/MKL/BLIS 1. Gradient and energy tolerances
are `1e-3` and `1e-7` Eh respectively. Original input iteration limits are
retained; no additional macro-iteration cap is set. Slurm requests a two-hour
allocation per pair and 64 GiB per process allocation.

Native wall time and peak process RSS are measured with GNU `time -v`.
HVP counts come from accepted outer-iteration NEO trace rows. On an incomplete
or failed run, these counts may omit work on the unfinished final macrostep;
they must not be presented as complete total action counts. Final energy and
explicit VBSCF convergence are parsed separately from Slurm exit status.

## Reporting

`summarize_preconditioner_suite.py RUN_ROOT` prints all 30 records as CSV.
Its `--markdown` output includes full energies and paired energy differences.
Pending/running entries retain missing final values, rather than being treated
as zero-time or unconverged completed runs. The raw logs, NEO traces, commands,
input/binary hashes, CPU/node and thread settings are retained under
`preconditioner-suite-1247682/` locally and in the remote run root.

A nominally converged pair with different final energies cannot be ranked
as solving the same optimization problem successfully. In particular,
FeCl2 must be checked against the established low-energy endpoint, not only
against the optimizer's convergence message. Small single-run wall-time
differences are not evidence of statistically significant speedups.

This selected suite can establish performance on these inputs; it cannot
establish universal superiority across all HAO/OEO wave functions.
