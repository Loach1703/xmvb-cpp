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
