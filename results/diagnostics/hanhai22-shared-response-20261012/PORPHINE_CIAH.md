# Fe porphine: current CIAH convergence runs

Submitted Hanhai22 Slurm array **1247745** on 2026-10-12:

| Task | State | Multiplicity | Active electrons/orbitals | PySCF final energy (Eh) |
|---|---|---:|---|---:|
| 1247745_0 | singlet | 1 | 8/11 | -2244.734377710447 |
| 1247745_1 | triplet | 3 | 8/11 | -2244.782671264980 |

These are new convergence runs, not another four-step diagnostic. The original
`itmax=1000` is retained; no extra outer-iteration limit is imposed. Slurm uses
the partition's 15-day resource limit, not a short chemistry timeout.

Both run exact libcint, Cartesian cc-pVDZ (469 basis functions), FULL structures,
full-AO OEO, Davidson and single-state coupled CIAH/NEO with the explicitly
selected **analytic Hessian-diagonal** preconditioner. Thresholds remain
`1e-3` gradient and `1e-7` Eh energy change. Each has 32 physical-core workers,
cores/close binding, BLAS one thread and a 192 GiB allocation. The existing
256 MiB exact-exchange workspace setting is retained from the earlier runs.
No additional convergence or endpoint-stability criterion is introduced.

## Source and inputs

The latest committed tree at submission is `a99a86ef`. Its production exact
path is unchanged from the tested `ab7e3c09` implementation; subsequent commits
only add benchmark scripts, analysis and documentation. Therefore this run
reuses the validated frozen binary rather than rebuilding identical code:

- Binary: remote `ciah-positive-build-1247629/xmvb-cpp.exe`.
- SHA256: `10430688c99848326f6c21967d4bf7e43c4a2ba1db1628b6e7a63b1e064e98a1`.
- Script validates this hash before execution.
- Unrelated uncommitted RI source/test changes are not in this exact-integral
  executable. They are preserved locally and are not silently committed.

Inputs are the previous corrected shared occupied-preserving ROHF frame,
with coefficients embedded in `$GUS`, not generated afresh. Remote and local
input hashes were checked before submission:

- Singlet: `e7d58d653a5d99f6acee301e0096e582dedd7a9f47c22b4d2c2ace013ab351fb`.
- Triplet: `79edd76ddc6902f84492592cfd70054a449b403ea150b4429ce5306e415ce8eb`.

Local input copies:
`porphine-converge-1247233/{singlet,triplet}/input.xmi`.
Remote input source: `inputs/{singlet,triplet}/input.xmi` in the run root.

## Reference and reporting

Existing PySCF 2.14.0 references are reused, not recomputed:
`results/studies/fe-porphine-current-20261012/pyscf-reference-284019/`.
They use the same saved initial orbital frame, Cartesian basis, 8/11 active
space, charge, spin target and no molecular symmetry constraint. Recorded
spin-square values are approximately 0 and 2 for singlet and triplet.

PySCF initial energies are -2244.704741852711 and -2244.611861070039 Eh.
Historical CPP initial energies differed by approximately 1.69 microhartree;
this difference must be reported rather than assuming exact initial agreement.
No additional fixed-point energy rerun is requested in this experiment.

Output root:
`/home/guqqgroup/xiatao/xmvb-runs/shared-response-20261012/porphine-ciah-1247745/`.
Each state saves the full output, NEO trace, native GNU-time wall/RSS, command,
binary/input hashes, CPU/node, build flags, thread settings and an affinity
probe. Later comparison must include initial/final energies, signed difference
from PySCF, outer steps, HVP count and native wall/RSS. Scheduler completion or
an optimizer convergence message alone does not establish energy agreement.
