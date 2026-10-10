# MnF2 multicore follow-up: shared spin-string links

## Finding and change

The previous dynamic AO-loop patch remains in place. This audit targets the
selected-direction phase that became slower with more threads. No MPI,
optimizer changes, relaxed tolerances, or molecule-specific thread thresholds
are introduced.

Instrumentation of the unchanged implementation (job 280811, Hanhai25/anode009)
isolates the cost:

| Threads | Selected direction / s | Connection construction / s | Connection builds |
| ---: | ---: | ---: | ---: |
| 8 | 0.124 | 0.014 | 1376 |
| 16 | 0.647 | 0.504 | 1376 |
| 32 | 2.595 | 2.225 | 1376 |

Every integral-direction sigma object rebuilt the same spin-string excitation
links, including fresh small OpenMP regions and allocation/hash lookup work.
These links depend only on the ordered occupation lists, not on the orbital
coefficients or integral values. Sharing them therefore changes no matrix
element, derivative, or optimization space.

DirectCiSigmaAction now owns immutable shared alpha/beta link tables. A second
constructor reuses those links and builds only the new numerical Hamiltonian
kernel. Both selected-direction and adjoint-direction consumers use it. Equal
spin spaces share the alpha table; distinct spin spaces retain separate tables.
Shared ownership preserves validity after the original sigma object dies.
Unused occupation-list copies and their accessor in OrthogonalDirectCiData
were removed. There are no new thread-private operator, integral, or carrier
copies, and existing bounded sigma workspaces are unchanged.

Across 669 HVPs, connection builds fall from 1376 to 38, removing 1338 small
parallel construction regions. Remaining builds occur when establishing the
fixed-point sigma action. They are not skipped using an accuracy threshold.

## Matched full-program comparison

Job 280816, anode009, exact Libcint/def2-TZVP MnF2 HAO, NEO/Davidson.
The unchanged ISCF=7 input is explicitly overridden by --optimizer-backend neo
in both runs, as in the preceding audit. Gradient tolerance 1e-3 and energy
tolerance 1e-7 hartree; BLAS threads=1, compact allocated-core affinity,
OMP_PLACES=cores, OMP_PROC_BIND=close, no wait-policy override.

Baseline is the immutable instrumented pre-sharing executable from job 280811.
Both versions have identical phase instrumentation. The first capture measures
the shared-link patch; final removal of unused occupation-list fields and
additional lifetime/equal-spin tests were subsequently verified in job 280817
and included in the reverse-order full runs below.

| Threads | Before / s | Shared links / s | RSS before / MiB | RSS after / MiB | Efficiency before | Efficiency after |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 161.436 | 159.693 | 383.8 | 383.1 | 100% | 100% |
| 8 | 26.701 | 26.963 | 446.5 | 428.1 | 75.6% | 74.0% |
| 16 | 22.101 | 17.302 | 469.5 | 449.4 | 45.7% | 57.7% |
| 32 | 19.411 | 14.600 | 517.4 | 497.8 | 26.0% | 34.2% |

Time is application end-to-end wall time; RSS is GNU time process peak.
Efficiency uses each version's measured single-thread time divided by its
parallel time and thread count. All eight runs converge in 14 iterations,
with 669 orbital HVPs. Energy range: -1348.893353119325 to
-1348.893353119304 hartree, a spread of about 2.1e-11 hartree.

The 16/32-thread wall reductions are 21.7%/24.8%. The 8-thread capture is
1.0% slower. Selected-direction time at 32 threads falls from 2.491 to
0.277 seconds; at 16 threads, from 0.687 to 0.119 seconds. Other AO-phase
times also improve despite their kernels being unchanged in this patch.
The whole-program gain must not be presented as a theoretical consequence
of the selected-phase savings alone: those additional effects are measured,
but their allocation/cache/runtime mechanism has not been isolated.

## CPU accounting versus useful parallel work

| Threads | Mean process CPU cores before | After | After / allocated cores |
| ---: | ---: | ---: | ---: |
| 8 | 7.65 | 7.65 | 95.7% |
| 16 | 15.19 | 14.83 | 92.7% |
| 32 | 28.21 | 28.80 | 90.0% |

These are GNU time user+system CPU divided by GNU elapsed time. They include
OpenMP spinning, so they are CPU busy/accounting indicators, not useful-core
counts. For example, 16-thread efficiency improves even though average
accounted cores decrease slightly, because removing redundant work and waits
is valuable even when it reduces CPU consumption.

The main AO loops retain roughly 97% loop wall-time occupancy/balance at
16/32 threads. This is a different metric: worker loop time divided by worker
capacity over the critical loop span, excluding the following barrier. It
includes memory stalls and scheduling delays, and is not whole-program
efficiency. Phase process-CPU counters are retained as diagnostics only;
short-window clock accounting can even yield CPU/wall ratios above the thread
count. Do not interpret those ratios as precise per-phase utilization.

## Reverse-order control and memory qualification

Job 280819 runs the final cleaned implementation first, then the unchanged
baseline for each core count, on anode033. Each matched pair has the same
allocated CPU set. This checks execution-order/background-load effects;
it is not a three-repeat timing campaign. No single-thread calibration was
repeated on this second node, so no efficiency is inferred for this capture.

| Threads | Baseline / s | Final / s | Baseline RSS / MiB | Final RSS / MiB |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 25.579 | 24.254 | 417.5 | 442.2 |
| 16 | 20.517 | 15.836 | 464.8 | 471.9 |
| 32 | 17.948 | 13.382 | 504.4 | 513.3 |

The large 16/32-thread gains persist (22.8%/25.4%), and the same connection
count reduction recurs. Eight-thread results are not a robust speedup claim.
RSS varies and increases in this control: removing duplicate logical data is
not a demonstrated reduction in process peak RSS. No large arrays were added
per thread; no universal RSS improvement is claimed. All six control runs
again converge in 14 iterations with 669 HVPs and numerically matching energy.

## PySCF reference and remaining targets

PySCF's [FCI contraction implementation](https://github.com/pyscf/pyscf/blob/master/pyscf/lib/mcscf/fci_contract.c)
shares compressed excitation links across workers and uses a parallel region
around tiled scatter/GEMM/reduction, with bounded per-worker tile buffers rather
than per-worker integral copies. This informs the separation of immutable
topology from numerical kernels here; its fixed block constants and orthogonal
orbital optimization rules are not copied.
Its [NumPy helper implementation](https://github.com/pyscf/pyscf/blob/master/pyscf/lib/np_helper/np_helper.c)
also illustrates first-touch initialization and contiguous SIMD loops. Those
are leads for future work, not changes implemented by this patch.

The 32-thread capture still spends about 3.15 seconds in NEO projected solves,
1.32 seconds in metric solves, and substantial time in exact AO contractions.
Further work should separate short serial/projection sections, repeated
parallel-region transitions, and shared-data streaming before adding threads
or changing workspace sizes.

Job 280820 provides one bounded full-MnF2 gprofng capture at 16 threads. It
resolves AVX-512 contract_rows, fused/transpose AO-H1E kernels, and OpenMP
barriers as remaining leads. Hardware and software perf probes are denied.
The collector again warns that its timer period changed: sample percentages
are not reliable quantitative costs. Neither profile time nor sampled barrier
percentages are used to claim native speedups or removable wall time.

## Tests and artifacts

- Job 280816: active-two-electron sparse, exact/RI AO-H1E, orthogonal direct-CI
  tests at 1/8/32 threads: nine passed executions.
- Job 280817: expanded direct-CI tests at 1/8/32 threads: three passed,
  including fresh-versus-shared sigma, owner lifetime, and equal-spin sharing.
- Job 280819: the same three suites at 1/8/32 threads: nine passed executions.
- Build flags remain -O3 -DNDEBUG -g1 -fno-omit-frame-pointer, Eigen BLAS on.
  Per-run hashes, CPU masks, resources and outputs are retained.

Remote root: /home/guqqgroup/taoxia/xmvb-runs/native-20261011/

Local raw evidence: results/studies/locality-20261011/; tables ci-sharing.csv,
ci-sharing-reverse.csv, and ci-sharing-phases.csv. Raw runs are not committed.
The profiler capture and executable remain on the cluster.

locality_results.py extracts converged native runs and leaves speedup/efficiency
blank if no single-thread calibration exists. phase_results.py writes phase
and loop-balance CSVs without summing overlapping scopes. mnf2-storage.sbatch
supports a reverse mode for paired order controls; mnf2-ci-verify.sbatch verifies
final source cleanup on compute nodes.
