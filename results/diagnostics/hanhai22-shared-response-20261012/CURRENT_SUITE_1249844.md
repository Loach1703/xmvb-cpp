# Current NEO multi-system validation, 2026-10-13

## Execution

Hanhai22 array **1249844**, 15 independent cases, 32 bound cores and 64 GiB
per process. Exact Libcint integrals, Davidson, analytic Hessian diagonal,
the unchanged 1e-3 gradient / 1e-7 Eh energy thresholds, and a 256 MiB
exact-exchange workspace budget. Original inputs are copied unchanged; no
30-macro diagnostic cap is added. Slurm uses the partition maximum time.
No old binary or molecular baseline is rerun.

Frozen restored build: **1249820**, all five unit executables passed.
Executable SHA256:
`9f51d9c85e5fabe3dcfd11108ff2075e86b60bb883b40b55d0922a8c70e2ab01`.
The active algorithm is the retained balancing CIAH implementation at
096d2419, not either rejected response-aligned/complement-completion trial.
Opt-in origin finite-difference auditing is disabled; scalar phase timers
do not add curvature or inverse actions. Local uncommitted RI changes are
not included in this frozen exact-integral build.

`summarize_current_suite.py` parses outputs, GNU time and NEO counters.
Historical comparisons use `preconditioner-suite-1247682` only after input
SHA256 equality and confirmed successful convergence. Older nodes, builds
and workspace budgets differ; wall differences are not matched speedups.
RSS is GNU process peak, not compiler memory or Slurm polling. On failures,
NEO TSV counters are accepted-step subtotals, not all attempted work.

## Preliminary snapshot

12 cases converged: F2, 241, 240, MnF2, FeCl2, 7975, 7963, YAMSAI,
NEGPIA, QOCGUM, 67970 and 10698. F2-OEO failed; CERRAS and LOFLEA remain
running at this snapshot. The generated `current-suite-1249844/summary.csv`
is explicitly partial, not a declaration that all jobs completed.

FeCl2: 12 macros / 248 HVPs, -2181.624464871811 Eh, 8.27 s, 491.53 MiB.
It is 1.1922e-8 Eh above the strict PySCF endpoint. Other converged results
agree with their historical energies to better than 7e-10 Eh, except FeCl2
whose old diagonal result was less accurately minimized. All completed
successful cases have fewer HVPs than the saved same-input comparison.
This does **not** mean all are faster: NEGPIA, 67970 and especially QOCGUM
have larger wall times, requiring attribution before improvement claims.

## Two issues found

**F2-OEO zero-macro exit.** The first macro is accepted. At macro 2's first
usable AH candidate, a nonlinear keyframe mismatch rejects the candidate.
The code subtracts that sole microstep, leaving zero accumulated displacement,
then exits with `ciah_zero_macro_step` without trying further AH enrichment.
Full gradient is still 0.00402497, so this is not convergence. The candidate
keyframe had lower energy, but the unchanged mismatch criterion rejected its
gradient model; it cannot simply be admitted because its energy decreased.
The next repair should address zero-displacement rejection/continuation,
not weaken the mismatch or convergence tolerance. The tested binary remains
frozen; this test does not silently introduce a repair.

**Selected-structure metric inverse cost.** QOCGUM converges in 3/37 versus
old 4/55, yet takes 52.32 s versus the historical 15.00 s. Its complete
structure norms cost **41.3163 s**, and their CG solve reaches **557**
iterations. Expansion HVP scope costs 1.8037 s and AH spectrum/image work
0.00445 s. These scopes can nest and must not be added into an invented
wall total. Full metric certification is an actual bottleneck here.
The partial CERRAS and LOFLEA captures likewise show hundreds of iterations
and over 100 s in structure norms. These are not completed runtime totals.
FeCl2's one-iteration metric behavior cannot be generalized to selected HAO
structure spaces. No new precision threshold or stopping condition is used
to disguise this cross-system weakness.
