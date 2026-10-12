# Porphine original-guess reference status

2026-10-13 observation on Hanhai22. All use the original supplied orbitals,
CAS(8,11), 469 Cartesian cc-pVDZ AOs, no ECP, and 32 allocated CPUs.
The earlier geometry-conversion hypothesis was disproved; do not revive it.

| Spin / method | Job | State | Macro | Energy / Eh | GNU wall | RSS |
|---|---|---|---:|---:|---:|---:|
| Triplet PySCF default AH | 1248772_3 | Converged | 35 | -2244.7908349995173 | 2667.29 s | 83.16 GiB |
| Triplet NEO | 1248745_2 | Running | 30 observed | -2244.7826442644 | Pending | Pending |
| Singlet PySCF default AH | 1248772_1 | Running | 30 observed | -2244.73272213814 | Pending | Pending |
| Singlet NEO | 1248745_0 | Running | 39 observed | -2244.7252037002 | Pending | Pending |

PySCF reports 639 JK and 117 micro iterations for the completed triplet;
those counts are not equivalent to NEO HVP counts. HVP count was not recorded
for this reference. Its result.json reports converged=true and spin square
2.0000000000275917; Slurm COMPLETED, exit 0. Its timer scope inside Python
is 2664.76 s, distinct from end-to-end GNU wall above.

The NEO jobs use the frozen d6f195c1 executable from build 1248744, before
the composed frame-component correction of 6d1f4853 and the later diagnostics.
They were not silently updated during execution. Neither NEO endpoint is
established yet, so the goal's numerical equivalence and efficiency comparison
remain unproven. No replacement porphine run was submitted this round.

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
results. The existing long jobs and the PySCF singlet are still running;
terminal energy equivalence remains unverified. Next, investigate whether
the already computed pure-orbital/pure-structure images should enrich the
actual AH trial space, rather than only its residual preconditioner.
