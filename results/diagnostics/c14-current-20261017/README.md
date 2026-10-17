# C14: current qualified implementation versus saved historical convergence

Hanhai22 job **1265059** completed successfully on **anode42**, 32 physical
cores, affinity **96-127**, with `ciah_projected_gradient_tolerance`.
Only the current implementation was submitted; neither GitHub nor the old
C14 native executable is rerun.

Production revision: **c45f9942** (numerical implementation d54c9095).
The rejected projected-channel experiments are absent. The execution mirror
passed the source manifest from qualification job **1265047**, which rebuilt
the restored implementation and passed 11/11 regressions. The benchmark job
then incrementally built and froze its own executable. Unrelated local RI
changes are not included.

## Input and native measurement

Original C14 anthracene geometry and guess, singlet CAS(14,14), STR=COV,
429 structures, HAO, exact Libcint/cc-pVDZ, ISCF=9/CIAH and Davidson.
Input ITMAX=200, convergence thresholds, coordinates and response policy
are unchanged. The original title mentioning four Kekule structures does
not override STR=COV. No profiler or extra macro-step limit is enabled.

- Executable SHA256: `ec71ee07f7e3677f4fed774cde6130c4d8af392544478275309deacec55bd057`.
- Input SHA256: `c9600592898c7cadf087d9cb3f77854a83eddcccc51962ac3c20bb0f393520b3`.
- Basis SHA256: `d53e7e02d55fb99cdf823d50cf13fec375e0cbe1220ff61816aee7ab38e38ec9`.

OpenMP uses 32 workers, close/cores binding, and one BLAS thread. Native
application wall time, GNU-time process wall/peak RSS, iteration output and
per-step HVP/response counters are retained separately. Build time is
not included in the application timer. GCC14.1/C++23 Release configuration
and source/binary/input hashes are retained by the batch script.

Remote artifact directory:
`/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015/c14-current-1265059`.

## Saved historical control, not rerun

Job **1264476**, anode40, 32 cores, affinity60-91:

| Quantity | Saved result |
|---|---:|
| Initial energy / Eh | -534.394411525272 |
| Final energy / Eh | -535.672671742972 |
| Accepted steps | 4 |
| Orbital HVP / structure responses | 46 / 46 |
| Per-step HVP | 5,10,28,3 |
| Application wall / s | 1190.925493 |
| Process wall / s | 1191.54 |
| GNU-time peak RSS / KiB | 6400232 |
| GNU-time peak RSS / GiB | 6.1037 |

The historical node was fully allocated at submission, so the new run used
anode42. Both nodes report AMD EPYC 9654 CPUs; the comparison still has a
node/placement difference.
Historical source and evidence: `../c14-ciah-20261016/metric-fix.md` and
`../c14-ciah-20261016/fix-1264476/`.

## Completed native comparison

| Quantity | Saved historical | Current retained |
|---|---:|---:|
| Initial energy / Eh | -534.394411525272 | -534.394411525272 |
| Final energy / Eh | -535.672671742972 | -535.672671742975 |
| Accepted steps | 4 | 4 |
| Orbital HVP / structure responses | 46 / 46 | 46 / 46 |
| Per-step HVP | 5,10,28,3 | 5,10,28,3 |
| Application wall / s | 1190.925493 | 306.172876 |
| SCF wall / s | 1184.541461 | 300.213401 |
| Process wall / s | 1191.54 | 306.62 |
| GNU-time peak RSS / KiB | 6400232 | 6078672 |
| GNU-time peak RSS / GiB | 6.1037 | 5.7971 |

The inputs are byte-identical. No metric-CG stagnation occurred. Application
time decreases 74.29% (3.89x observed acceleration); peak RSS decreases 5.02%.
The unchanged per-step action counts distinguish lower internal cost from
fewer HVPs or accepted steps. These are single native historical comparisons,
not a same-allocation A/B or an attribution of every saved second to one kernel.

`evidence.tar.gz` contains the completed current run (without its executable)
and the historical native output, input, trace and resource records. It can be
read directly by `benchmarks/history/summarize_c12_c14.py`, without extraction
or running an old executable. The combined report and CSV are in
`benchmarks/history/c12_c14_comparison.md` and `.csv`.
