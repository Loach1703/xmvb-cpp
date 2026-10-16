# COV: eliminate repeated exterior maps in the selected metric solve

## Bottleneck and necessary algebra

This round changes neither the CIAH policy nor its thresholds. The selected
metric and CI transport are necessary; repeating the same exterior maps in
every CG action at a fixed chart is not.

Current-only attribution job 1263201 found six selected metric solves with
1,684 operator applications, taking 15.914474 s. Their action time was
15.909127 s. Overall left exterior transforms were called 112,400 times
(17.976488 s); the six orbital HVPs took 2.038697 s. These are explicit
nested phase clocks, not additive sampled-profiler percentages.

Additional panel attribution 1263214 found 210 ordinary product H/M actions,
3,150 input panels, 3,075 output projections, and **zero input/output
evictions**. Increasing tile residency would not address the dominant loss.

For singleton beta columns, including the complement matching in covalent
spaces, the raw scalar coefficient can be moved outside both exterior maps.
The coefficient-free alpha/beta basis images are shared across a RHS batch.
Their sampled product is contracted with all RHSs by GEMM and finally with
the sparse structure map. This is exact for arbitrary selected structures
on that support; there is no closure assumption or amplitude truncation.

The fixed selected metric is then projected once when its selected-space
representation fits the existing workspace. Subsequent CG actions apply
that same metric without returning to the spin-string carrier. This changes
only the fixed linear operator's storage/contraction order, not the physical
representation or its admission to orthogonal CI. The existing CG, initial
vector, preconditioner and stopping conditions are retained. The metric is
rebuilt whenever its orbital chart changes; its derivatives still use the
existing analytic overlap adjoints.

Section 10.5 of `article/source_column_hamiltonian.tex` records the equations,
work and storage bounds. The 25-page PDF was compiled and visually checked.

## Memory

The selected-metric workspace remains 8 MiB: half for construction, one
quarter each for the selected identity and result. No full string-pair Gram
or transformed spin-product array is introduced. The 132-structure COV metric
occupies 139,392 bytes (0.133 MiB). Source/Gram panels and RHS product blocks
are jointly bounded. Larger selected metrics retain their exact tiled
operator. Unneeded exterior objects are released after the small selected
metric has been projected. No per-thread carrier or integral array is added.

This storage choice does not claim to solve the large-selected-metric case;
it removes a demonstrated repetition for the present COV workload without
making large spaces allocate quadratic string-pair data.

## Native end-to-end one-step results

LOFLEA CAS(12,12), cc-pVDZ, HAO, libcint exact integrals, NEO/CIAH,
structure-space Davidson, analytic Hessian diagonal, 1e-3 gradient and
1e-7 Eh energy tolerances. Eight actual OpenMP cores within one LLC,
single-threaded BLAS; 16 CPUs allocated and 48 GiB requested. COV has 132
structures and FULL 226,512. GNU-time process wall includes final output.

| Space/version/job | Node | SCF (s) | Process wall (s) | RSS (GiB) | HVP/response | Final CI actions |
|---|---|---:|---:|---:|---:|---:|
| COV saved remote / 1262622 | anode50 | 53.464025 | 58.06 | 2.612 | 6/6 | 99 |
| COV previous current / 1263149 | anode04 | 39.691606 | 43.61 | 2.518 | 6/6 | 99 |
| COV final / 1263242 | anode04 | 30.096030 | 34.55 | 2.531 | 6/6 | 99 |
| FULL saved remote / 1262623 | anode50 | 74.965957 | 93.77 | 3.373 | 8/8 | 39 |
| FULL previous current / 1263147 | anode04 | 54.074424 | 75.75 | 3.350 | 8/8 | 27 |
| FULL regression / 1263227 | anode04 | 54.225092 | 73.14 | 3.242 | 8/8 | 27 |

COV process wall falls **20.8% against the previous same-node current** and
40.5% against the saved remote. The remote comparison crosses nodes with
the same EPYC 9654 CPU model; no remote run was repeated. Different core
groups and concurrent node load remain timing variables. Previous-current
COV used CPUs 176--183; final affinity is retained in its JSONL artifact.
No repeated-run average is claimed.

The first qualification binary (before releasing the now-unused transform
objects and removing one overwritten primal expression) used 36.28 s and
2.495 GiB in job 1263225. It is an implementation-stage qualification, not
another independent sample to average with the final binary. RSS is roughly
unchanged at 2.5 GiB, not a demonstrated end-to-end memory reduction.

FULL does not enter the new small selected-metric representation. Its
unchanged SCF time/actions and energy are a regression check; its process-wall
and RSS variation are not attributed to this COV optimization.

## Numerical checks

| Space | Initial energy (Eh) | After one step (Eh) | Difference from previous current (Eh) |
|---|---:|---:|---:|
| COV | -398.533053417516 | -399.074297727156 | +3.2909e-8 |
| FULL | -398.880685140388 | -399.322366416763 | +4e-12 |

COV differs from the saved remote after one step by +2.9935e-8 Eh. All runs
are **one accepted optimization step, not converged energies**. Exit 2 with
`max_iterations` is expected. COV retains 10 microsteps and three keyframes;
FULL retains 13 and three. HVP and final CI counts are unchanged. Reordered
metric arithmetic changes finite-precision CG iterates/counts, but no
residual tolerance or optimizer test is relaxed.
Final physical CI residuals remain below the unchanged 1e-7 Eh budget:
7.65005318e-8 Eh for COV and 6.34373532e-8 Eh for FULL.

Final build 1263235 passed 12/12 CTest checks on Hanhai22, including independent
occupied-minor cross-Gram references, chart directional differences, frozen
prepared metrics, CI transport and keyframes. New singleton-column cases
cover unequal and equal alpha/beta dimensions and arbitrary nonclosed
three-column selections, including tiny-workspace comparisons.

Additional CAS(12,12) verification 1263255 compares the prepared selected
metric with independent coefficient-weighted, one-RHS-at-a-time factor
actions for 1, 3, 20 and 132 selected covalent structures. Relative errors
are respectively 1.96e-16, 4.64e-16, 3.48e-16 and 5.98e-16. Existing
duality and curvature finite-difference checks also pass. This correctness
job ran on anode31 after anode04 filled; it is not a molecular timing sample.

## Attribution and remaining cost

| Explicit phase | Before / 1263201 (s) | After / 1263226 (s) |
|---|---:|---:|
| Six selected metric solves, including projection setup | 15.914474 | 0.072960 |
| Left exterior maps | 17.976488 | 5.425954 |
| Two keyframes | 14.645170 | 2.715819 |
| CIAH setup | 3.168465 | 0.112240 |
| Orbital HVPs | 2.038697 | 2.067272 |

The after-attribution binary precedes the two final cleanup changes described
above; its mathematics and thresholds are identical. Instrumented timings
are not used for production speed claims. The final optional attribution
job 1263256 never ran and was cancelled after anode04 filled; no queued
benchmark is left behind.

The remaining bottleneck is the initial and trial CI solution: COV still
needs 99 final H/S actions, whereas FULL uses 27. Fewer selected coefficients
do not eliminate the orthogonal intermediates of a nonclosed space. Neither
that necessary support nor the selected projection may be dropped. The next
useful analysis is the conditioning/work of this selected CI solve, rather
than HVP budgets, larger tiles, or additional convergence certifications.

## Evidence

`cov_selected_metric_scf.csv` combines saved remote/previous-current records
with the final native results. `cov_selected_metric_phases.csv` holds the
phase extraction. The evidence archive retains input/binary hashes, affinity,
raw phase logs, traces, GNU-time and Slurm data, final build/CTest records,
CAS12 parity outputs and the scoped source patch. Full native FULL stdout
remains on Hanhai22; its extracted summary and SHA-256 are included locally.
Unrelated local RI edits were not included in this change or transferred.
