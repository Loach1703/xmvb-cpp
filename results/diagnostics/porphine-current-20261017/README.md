# Current iron-porphine timing runs

Current production code: `c491230b`. Single-state FULL OEO, exact Libcint,
CAS(8e,11o), cc-pVDZ, original per-spin XMI guesses, CIAH/Davidson and
analytic Hessian diagonal. The gradient and energy tolerances remain
1e-3 and 1e-7 Eh. No shortened SCF run was introduced.

## Migration to Hanhai25, 2026-10-17

Hanhai22 array `1265637_[0-2]` stayed pending on scheduling priority. It was
held during the replacement environment check, then cancelled without having
started (confirmed CANCELLED and elapsed 00:00:00). Old results are preserved.

The exact qualified executable, Libcint/OpenBLAS libraries, basis and inputs
were copied to Hanhai25; no old code was executed and no rebuild was needed.
The executable SHA256 remains
`7f1aec8c5a36661027427f92c667e83d3e0c5ed7582aba0098188e50b146033e`.
The original build is GCC 14.1.0/C++23; the Hanhai25 runtime loads GCC 14.3.0.
Qualification `291576` completed successfully on anode006: dependency
resolution, product-tile unit regression and converged 32-core 241 smoke test.

| Spin | Hanhai25 task | Node at start | Migration check |
|---|---|---|---|
| Singlet | 291577_0 | anode006 | RUNNING |
| Triplet | 291577_1 | anode022 | RUNNING |
| Quintet | 291577_2 | anode029 | RUNNING |

Each is one process with 32 physical cores, cores/close binding and BLAS1.
The 192 GiB Slurm request is an allocation limit, not measured peak RSS.
GNU time and Slurm accounting will supply actual resource measurements.
The qualification CPU is AMD EPYC 9755, whereas the earlier Hanhai22 runs
used EPYC 9654. Cross-cluster wall ratios must therefore retain that distinction.

Raw results:
`/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/run-291577/`.
Original submission metadata is in `submission.json`; current migration and
job mapping are in `migration-hanhai25.json`. The two batch scripts preserve
the original and replacement environments. Qualification evidence and the
new bundle manifest are archived locally. Historical code and other software
are not rerun.
