# Restored CIAH strategy: molecular regression

Submitted on Hanhai22 on 2026-10-14. All jobs use the source at master
commit `39201b81`, frozen by successful build/test job `1255480`.
Binary SHA256: `a5313100b12fd57e7d876e0d21c9333df5cc9d1a4d4dfca126c25b3f5da9b602`.
Uncommitted local RI modifications are excluded.

| Job | System |
|---|---|
| 1255488_0 | Fe-porphine singlet (previously submitted) |
| 1255492_0 | F2-OEO |
| 1255492_1 | FeCl2 |
| 1255492_2 | MnF2 |
| 1255492_3 | 241 |
| 1255493_2 | Pyrrole dimer OEO, two-state average |
| 1255493_3 | Pyrrole dimer BLO, two-state average |

Use 32 cores per task, exact integrals, Davidson, analytic Hessian-diagonal
preconditioning, original inputs, gradient tolerance 1e-3 and energy tolerance
1e-7 Eh. No added short-run iteration cap. The pyrrole array uses the existing
`pyrrole-recheck-20261014/ciah.sbatch` with `BUILD_JOB=1255480` and indices 2-3.
No historical version is rerun. Completion and numerical equivalence must be
checked before reporting performance. Record final energy, accepted macros,
orbital HVP, GNU wall time, peak RSS, node and binary/input hashes.
