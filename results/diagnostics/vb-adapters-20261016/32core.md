# Current 240/241: 32-core native follow-up

Both jobs converged on Hanhai22 `anode38`, launched concurrently as array
`1264434_[0-1]`. No historical executable was run and no production code,
optimizer parameter, convergence threshold or input was modified for this test.

## Matched current runs

The same frozen C++23 Release executable used for the previous eight-core runs
was reused, SHA256
`27966b5ea621c92392457a56ecd6e25d87167e3f31c56039333b64d963d70b57`.
It represents committed production sources in `48c145ba`, with the earlier
448-file source-manifest qualification retained in `report.md` and
`evidence.tar.gz`. Pre-existing uncommitted RI worker-budget changes were not
included in this executable. No rebuild was needed.

Both inputs are byte-identical to the eight-core regression:

- `240_vb.xmi`: `61d1035fd059c52a7361d6ed8f49f00a125286c644c155127bf6f8ee865eba48`.
- `241_VBSCF.xmi`: `5dbd9c72c4ab5b8dec9b1bda5ef8c8d5d64741cc181f25cc8604e2464f9f07c8`.

Explicit CIAH/Davidson, exact Libcint integrals, analytic Hessian-diagonal
preconditioning, projected-gradient tolerance `1e-3`, energy-change tolerance
`1e-7 Eh`; no optimizer step cap added. OpenMP workers 32, cores/close affinity,
maximum active levels 1, BLAS threads 1. The actual application-step CPU masks
are `96-127` for 240 and `0-31` for 241, on a machine with one thread per core.
No optional phase profiler was enabled in either current run.

| Case | Threads | Steps | HVP | Initial energy / Eh | Final energy / Eh | SCF / s | Application total / s | Process / s | Peak RSS / MiB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 240 | 8 | 4 | 30 | -343.075645280883 | -343.446302545197 | 8.088943 | 10.963945 | 11.15 | 1528.16 |
| 240 | 32 | 4 | 30 | -343.075645280882 | -343.446302545196 | 2.856646 | 4.354682 | 4.53 | 1569.29 |
| 241 | 8 | 4 | 44 | -230.474471805007 | -230.720590533578 | 3.103752 | 4.170389 | 4.32 | 904.10 |
| 241 | 32 | 4 | 44 | -230.474471805007 | -230.720590533578 | 1.189183 | 1.916378 | 2.05 | 853.08 |

Both terminate with `ciah_projected_gradient_tolerance`, application and Slurm
exit status zero. HVP counts were summed from the per-macro `ciah.tsv` rows;
they do not change between current 8- and 32-core runs. Native process CPU
times are 99.85 and 43.88 seconds, approximately 22.0 and 21.4 utilized cores.
These are resource-accounted averages, not hotspot or thread-balance evidence.
RSS is GNU time's process peak; Slurm sampled RSS is uninformative for these
short jobs. Application total and external process wall time are separate.

## Saved historical comparisons

- **240, recent equal-energy Hanhai22 baseline:** saved job `1262045`, 32 cores,
  four steps, 30 HVP, SCF `3.082727 s`, application total `4.536277 s`, peak RSS
  `1570.26 MiB`. Current SCF is 7.33% lower, application total 4.00% lower;
  final energy differs by approximately `6e-12 Eh`. This does not show the large
  wall-time regression suggested by comparing the current eight-core run to
  the historical 32-core run. Hosts differ (`anode21` versus `anode38`), and
  these are single shared-node measurements, not statistical speedup claims.
  Source: `../orthogonal-tiles-20261015/full32_comparison.csv`.
- **241, earlier same-input history:** `../../history/neo_optimized_20261011.csv`,
  row `241_VBSCF`, records 27 HVP, SCF `1.374423 s`, total `1.822427 s`,
  RSS `918.63 MiB`, final energy `-230.720590367937 Eh` and the same input hash.
  Current SCF is lower, total is 5.16% higher, and 44 HVP remains more than 27.
  The historical endpoint is about `1.66e-7 Eh` higher and the solver revision
  differs; this is historical context, not a matched-refactor performance
  attribution. The approximately one-second `241_vb` records use a different
  initial guess and must not be substituted for the `241_VBSCF` input.

This follow-up separates the prior core-count mismatch from remaining
algorithmic work. It does not identify a particular contraction as a hotspot
or prove that the recent adapter cleanup caused the historical HVP difference.

`32core.csv` contains the current runs and the recent 240 baseline.
`32core-evidence.tar.gz` retains inputs, raw output, per-macro traces,
hash checks, affinity, environment, CPU/module records, GNU time and Slurm
accounting. It contains no executable or copied basis library.
