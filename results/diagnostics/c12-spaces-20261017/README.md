# C12 ionic, full and arbitrary selected-structure convergence

2026-10-17. Current implementation only; no historical executable was rerun.
These tests extend the COV-only evidence and do not change numerical code.

## Frozen setup

- C12 is LOFLEA CAS(12,12), HAO, 150 basis functions, exact Libcint/cc-pVDZ,
  singlet, one state. The original deck is
  `testdata/vbscf/csd_LOFLEA_vb.xmi`.
- Original geometry, strictly localized orbital supports, initial orbitals,
  Davidson and convergence criteria are retained. Only ISCF=9/CIAH and the
  requested structure space are selected. ITMAX remains the original 2000;
  no shortened optimization or application timeout is introduced.
- Source provenance is master **03d469cc**, numerical implementation
  **d54c9095**, qualified by **1265047** (11/11 tests). This reuses the frozen
  executable from converged C14 job **1265059**, not a new untested build.
  Executable SHA256:
  `ec71ee07f7e3677f4fed774cde6130c4d8af392544478275309deacec55bd057`.
- Hanhai22, **anode42**, 32 physical cores per task, one process,
  OMP_PLACES=cores / OMP_PROC_BIND=close, BLAS=1. Five array tasks were
  submitted together, not in a serial loop. Multiple tasks share this node.
  Each has a 64-GiB Slurm allocation; C++ workspace budgets are unchanged.
- No large numerical arrays are duplicated by this runner. All tasks execute
  one shared immutable executable. Inputs, basis and script hashes are saved.

Preparation **1265092** completed successfully; native array **1265093** runs
from `c12-spaces-bundle-1265092` in
`/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015` on Hanhai22.
The first preparation job 1265090 failed on a missing Python helper before
any molecular calculation; its dependent array 1265091 was cancelled. The
helper was included before the successful preparation and native submission.

## Structure spaces

| Array index | Label | Definition | Structures |
|---|---|---|---:|
| 0 | ion01 | ION(0:1) | 5676 |
| 1 | ion02 | ION(0:2) | 47256 |
| 2 | full | FULL | 226512 |
| 3 | mixed3 | Explicit covalent + first/second ionic structures | 3 |
| 4 | mixed20 | Explicit mixed covalent/ionic structures, orders 0–3 | 20 |

The explicit sets are fixed before calculation (seed 20261017), nested and
independent Rumer singlet structures within each occupation pattern. They are
not assumed closed under orbital transformations. Their exact `$STR` lines,
input SHA256 and generated-space counts are in `bundle-1265092/cases.json`.
All native runs report the expected structure count.

## Results snapshot

| Space | Status | Steps | HVP | Final energy / Eh | Application total / s | Process total / s | Peak RSS / GiB |
|---|---|---:|---:|---:|---:|---:|---:|
| ION(0:1) | Running | — | — | — | — | — | — |
| ION(0:2) | Running | — | — | — | — | — | — |
| FULL | Converged | 4 | 40 | -399.417495587054 | 97.036986 | 107.31 | 3.568096 |
| mixed3 | Converged | 4 | 38 | -399.050820856824 | 14.640989 | 14.95 | 2.662525 |
| mixed20 | Converged | 4 | 38 | -399.095231633359 | 20.185947 | 20.38 | 2.663704 |

Completed jobs have exit code zero and
`Termination reason: ciah_projected_gradient_tolerance`. Structure-response
action counts equal orbital-HVP counts. Initial energies, raw counts and
native timings are parsed into `summary.csv`, not transcribed from plots.
Running rows contain only the already accepted steps/actions, not a prediction
of total convergence work. The ION tasks remain submitted until completion.

Application total is the same native timer used in the historical C12/C14
report. GNU-time process total also includes remaining output/teardown work;
FULL has approximately 158 MiB of output at the original IPRINT=3. RSS is
GNU-time process peak, not Slurm's sampled maximum or allocated memory.

Different structure spaces have different variational energies. The completed
arbitrary mixed subsets demonstrate end-to-end applicability beyond COV, not
a measured speedup against a same-space historical run. ION is substantially
slower than FULL so far despite fewer structures; this record does not assign
a cause without phase timing or alter the solver to hide it.

## Existing historical comparisons

- Current COV (132 structures), job 1264965: 4 steps / 36 HVP,
  -399.106646267332 Eh, 17.764194 s, 2.656708 GiB.
  Historical same-space 32-core job 1264541: 4 / 36, same energy,
  44.996525 s, 2.597919 GiB. See
  `benchmarks/history/c12_c14_comparison.md`.
- Frozen GitHub revision 94693457 has a saved FULL record, but it is
  **8 cores and one accepted step**, not a converged 32-core control:
  8 HVP, initial -398.880685140388 Eh, endpoint -399.322429054264 Eh,
  application 80.973591 s, process 93.77 s, peak 3536844 KiB.
  Saved source: `results/diagnostics/orthogonal-tiles-20261015/github_scf_step.csv`.
  Do not report a matched full-run speedup from this pair.
- No saved matching C12 ION or mixed-subset convergence baseline was found.
  Existing 240/CAS(8,8) ION results are not controls for this CAS(12,12) task.
  These historical comparisons remain missing rather than rerunning old code.

## Artifacts and refresh

`evidence.tar.gz` contains the frozen input bundle (without executable) and
native result snapshots, including complete logs for finished tasks. The raw
authoritative live directory is
`/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015/c12-spaces-1265093`.
The archive is a snapshot, not an automatic live update.

Refresh only the result artifacts, then regenerate the CSV with
`python3 tests/mcscf/c12_spaces.py summarize PATH_TO_run-1265093`.
The preparation and submission procedures are implemented by
`tests/mcscf/c12_spaces.py` and `tests/mcscf/c12_spaces.sbatch`; they reuse a
qualified frozen binary rather than building or executing on the login node.
