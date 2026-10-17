# C14 CIAH convergence run

Submitted Hanhai22 job **1264455**, 32 physical cores, exact Libcint/cc-pVDZ,
HAO, singlet CAS(14,14), `STR=COV`, Davidson. The attachment's title mentioning
four Kekule structures does not override its actual `STR=COV` control.

`C14.xmi` differs from the supplied attachment only in `iscf=7` -> `iscf=9`.
Original SHA256: `7c6d6edf039007116b8944b52df7e8aacd961debcb94aa0638046282b5165efa`.
Run SHA256: `c9600592898c7cadf087d9cb3f77854a83eddcccc51962ac3c20bb0f393520b3`.
The original `ITMAX=200`, initial orbitals, geometry, basis and other controls
are unchanged; no additional optimizer cap or threshold override was added.

Executable SHA256:
`27966b5ea621c92392457a56ecd6e25d87167e3f31c56039333b64d963d70b57`.
This is the same qualified current-production binary as the 240/241 follow-up,
production-source revision `48c145ba`; no rebuild or historical run was needed.
Uncommitted local RI budget changes are not included.

Run directory:
`/home/guqqgroup/xiatao/xmvb-runs/c14-ciah-20261016/run-1264455`.
At launch, the job is running on `anode65`; convergence is not yet established.
Native output, per-macro CIAH counters, GNU time resource data, binary/input
hashes and actual core affinity are retained there. The standard performance
audit's Slurm/native-resource procedure is used on user-requested Hanhai22.
No optional profiler is enabled. OpenMP uses 32 workers and BLAS one thread.

The original run failed before accepting a step, at a selected-metric CG solve.
See [metric-fix.md](metric-fix.md) for the measured diagnosis, unchanged-accuracy
factor-reuse repair, passing regressions and replacement full run **1264476**.
