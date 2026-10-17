# C12 and C14: retained exact-integral performance

Updated 2026-10-17. The current fast numerical implementation is **d54c9095**;
the later c45f9942 qualification changes no numerical algorithm. The rejected
projected-channel experiments are absent from production. Restored-state
qualification **1265047** passed **11/11 regressions**. No old executable or
GitHub control was rerun for this report.

## Complete convergence, 32 physical cores

All four records use Hanhai22, exact Libcint/cc-pVDZ, HAO, a singlet covalent
VB structure space, ISCF=9/CIAH and structure-basis Davidson. OpenMP uses
close/cores placement and BLAS uses one thread. Initial orbitals, input files
and convergence criteria are unchanged within each molecular comparison:
projected-gradient infinity norm 0.001 and energy change 1e-7 Eh.

"C12" denotes the existing **LOFLEA CAS(12,12)** workflow, not a separate
carbon-C12 molecular geometry: 132 structures, 924 spin strings per spin,
150 basis functions. C14 is anthracene CAS(14,14), STR=COV, 429 structures.
Its input title mentioning four Kekule structures does not override STR=COV.

| System | Version / job | Steps | HVP | Final energy / Eh | Total / s | Peak RSS / GiB |
|---|---|---:|---:|---:|---:|---:|
| C12 / LOFLEA | Historical cd958ad7 / 1264541 | 4 | 36 | -399.106646267332 | 44.996525 | 2.5979 |
| C12 / LOFLEA | Retained d54c9095 / 1264965 | 4 | 36 | -399.106646267332 | **17.764194** | 2.6567 |
| C14 | Historical cd958ad7 / 1264476 | 4 | 46 | -535.672671742972 | 1190.925493 | 6.1037 |
| C14 | Retained c45f9942 / 1265059 | 4 | 46 | -535.672671742975 | **306.172876** | **5.7971** |

- C12: **2.53x** observed acceleration, **60.52%** less total time. Peak RSS
  increases by 60.20 MiB (2.26%); this is not a measured memory improvement.
- C14: **3.89x** observed acceleration, **74.29%** less total time. Peak RSS
  decreases by 314.02 MiB (5.02%). No metric-CG stagnation occurred.
- Per-step HVP counts are unchanged: C12 **6,17,10,3**; C14 **5,10,28,3**.
  Structure-response counts equal orbital-HVP counts in all four records.
- Initial energies are identical: C12 -398.533053417516 Eh; C14
  -534.394411525272 Eh. Final energies agree at numerical precision.

Total time is the native application's end-to-end timer, excluding build and
tests. GNU-time process times, which also include trailing process/output work,
are separately retained in the CSV. Peak RSS is GNU-time process peak, not
Slurm job memory, live heap, or the sum of thread-local allocations.

Historical records ran on **anode40**; current records ran on **anode42**.
Both report AMD EPYC 9654 CPUs. These are single historical comparisons across
nodes/placements, not repeat averages or a same-allocation timing experiment.
The action counts and energies support a lower internal cost, not a claim that
the measured speedup comes entirely from one function or parallelism alone.

## GitHub baseline is a separate timing scope

The frozen GitHub baseline is **946934570c096ba5072923effa7fb063736731f1**.
Its saved C12/COV record **1262622** used **8 cores** and only **one accepted
step**, not full convergence: 6 HVPs, application time **57.809377 s**, process
time **58.06 s**, peak RSS **2738996 KiB**. The endpoint is
-399.074297757091 Eh, not the converged energy in the main table.

Do not label the 44.996525 s historical convergence record as this GitHub
baseline, or divide the 8-core one-step time by a 32-core converged time and
call it a matched full-run speedup. Its original data remain in
`results/diagnostics/orthogonal-tiles-20261015/github_scf_step.csv`.

## Reproduce the table from saved records

Run only the Python parser, not any electronic-structure calculation:

```bash
python3 benchmarks/history/summarize_c12_c14.py
```

It reads the committed evidence archives directly, checks successful 32-core
convergence and identical input hashes/initial energies, and emits the full
precision CSV. `c12_c14_comparison.csv` is that output, with archive provenance,
job IDs, nodes, input hashes, native/process timings and action counts.

Evidence archives:

- `results/diagnostics/c12-tiles-20261017/evidence.tar.gz`: historical C12.
- `results/diagnostics/c12-gemm-orbit-20261017/evidence.tar.gz`: retained C12.
- `results/diagnostics/c14-current-20261017/evidence.tar.gz`: historical and current C14.

The detailed C14 report is
`results/diagnostics/c14-current-20261017/README.md`. The rejected C12 numerical
experiments and their explicitly unsuccessful performance results are archived
under `results/diagnostics/c12-projected-channels-20261017/`; none is a compiled
fallback or selectable production path.
