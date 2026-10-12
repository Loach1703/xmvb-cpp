# Porphine original-guess reference status

2026-10-13 observation on Hanhai22. All use the original supplied orbitals,
CAS(8,11), 469 Cartesian cc-pVDZ AOs, no ECP, and 32 allocated CPUs.
The earlier geometry-conversion hypothesis was disproved; do not revive it.

| Spin / method | Job | State | Macro | Energy / Eh | GNU wall | RSS |
|---|---|---|---:|---:|---:|---:|
| Triplet PySCF default AH | 1248772_3 | Converged | 35 | -2244.7908349995173 | 2667.29 s | 83.16 GiB |
| Triplet NEO | 1248745_2 | Energy rejection, not converged | 50 accepted | -2244.790272699886 | 5177 s | 38.63 GiB |
| Singlet PySCF default AH | 1248772_1 | Running | 62 observed | -2244.74099511752 | Pending | Pending |
| Singlet NEO | 1248745_0 | Energy rejection, not converged | 64 accepted | -2244.739920658597 | 5313 s | 38.62 GiB |

PySCF reports 639 JK and 117 micro iterations for the completed triplet;
those counts are not equivalent to NEO HVP counts. HVP count was not recorded
for this reference. Its result.json reports converged=true and spin square
2.0000000000275917; Slurm COMPLETED, exit 0. Its timer scope inside Python
is 2664.76 s, distinct from end-to-end GNU wall above.

The NEO jobs use the frozen d6f195c1 executable from build 1248744, before
the composed frame-component correction of 6d1f4853 and the later diagnostics.
They were not silently updated during execution. Both NEO jobs have now
terminated with `ciah_macro_energy_rejected`, Slurm exit 2. They are not
converged endpoints. Triplet differs from converged PySCF by +0.56230 mEh.
The goal's numerical equivalence and efficiency comparison remain unproven.
No replacement porphine run was submitted this round.

Triplet PySCF result, trajectory, resource record, node and hashes have been
copied to `porphine-rotation-1248772/triplet/pyscf-ah/` beside this note.

## Certified coarse-preconditioner four-macro test

1248959_0/1 completed the deliberately bounded four-macro tests (exit 2 is
expected, not a convergence claim). They use the frozen binary from 1248947,
the same original input orbitals and unchanged 1e-3/1e-7 tolerances.
No old executable was rerun. The old prefix is extracted automatically
from the ongoing original-guess 1248745 trajectories.

| Spin | Old prefix HVP | New HVP | Old macro-4 energy / Eh | New macro-4 energy / Eh | New end-to-end program time / s | New RSS / MiB |
|---|---:|---:|---:|---:|---:|---:|
| Singlet | 64 | 50 | -2244.7230503763 | -2244.7182260015 | 524.401 | 39546.1 |
| Triplet | 44 | 41 | -2244.7551098961 | -2244.7615530749 | 606.174 | 39559.4 |

This is mixed evidence: the singlet has fewer HVPs but is higher in energy
after four macros; the triplet has fewer HVPs and is lower. The old frozen
binary predates both the composed frame-step correction and the new coarse
inverse, so this is **not** an isolated preconditioner ablation. Matched
four-macro wall/RSS data are unavailable for that old running prefix. Do not
invent them or claim a paired timing speedup. GNU elapsed additionally
includes work outside the reported program end-to-end timer; preserve its
raw resource record separately.

No replacement unbounded NEO job is submitted based on these four-macro
results. Both old NEO jobs are now terminal; only PySCF singlet remains running.
Terminal energy equivalence remains unverified. The later split-Ritz experiment
did not improve FeCl2 and was removed; see `FECL2_SPLIT_AH.md`.

## Terminal rejection evidence

Singlet attempted macro 65 has an accepted keyframe at -2244.739935866974 Eh,
below source -2244.739920658597 Eh. Later frozen-curvature microsteps yield an
uphill full trial and the whole-step line search fails. Triplet attempted
macro 51 has its accepted keyframe at -2244.790247284049 Eh, already above
source -2244.790272699886 Eh. Final scaled directions in both cases have
positive origin derivatives, so shrinking toward zero cannot create descent.
This motivates examining nonlinear commitment and keyframe backtracking,
not invalidating an origin HVP at its own origin. Original-guess accepted HVP
totals are 667 singlet / 528 triplet; the final rejected macro's actions are
not included in those trajectory rows.
