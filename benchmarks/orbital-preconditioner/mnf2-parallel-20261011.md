# MnF2 native multicore audit, 2026-10-11

## Scope and result

Hanhai25 Slurm compute node anode009; one process, exact integrals,
NEO/Davidson, unchanged MnF2 input and convergence tolerances (gradient 1e-3,
energy 1e-7 hartree). No MPI and no RI changes. All matched runs converge
in 14 iterations with 669 orbital HVPs. Final energies lie between
-1348.893353119321 and -1348.893353119307 hartree.

The reproducible 16-versus-8-thread slowdown is primarily a static-loop
long-tail problem. Equal integral/CSR-edge counts do not imply equal elapsed
worker cost. Direct worker timers show a roughly 2.4-fold min/max duration
spread at 16 threads, with negligible start skew. Dynamically assigning
smaller shared pieces removes this slowdown. The precise hardware cause of
the unequal worker costs is not established: cache/frequency explanations
remain hypotheses, not measured conclusions.

## Matched final comparison

Job 280804 compares CSR-storage-ordered static loops (snapshot from job
280798) with the same implementation plus dynamic balancing. Both executables
include the same opt-in worker timing. Single captures, not repeated medians.
Each pair uses the same allocated CPU set and compact affinity policy.

| Threads | Static time / s | Dynamic time / s | Static RSS / MiB | Dynamic RSS / MiB | Static efficiency | Dynamic efficiency |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 158.902 | 158.587 | 382.9 | 383.0 | 100% | 100% |
| 8 | 26.346 | 27.202 | 423.2 | 422.4 | 75.4% | 72.9% |
| 16 | 31.275 | 22.490 | 459.6 | 472.3 | 31.8% | 44.1% |
| 32 | 22.413 | 19.725 | 504.9 | 508.0 | 22.2% | 25.1% |

Time is application end-to-end wall time. RSS is GNU time maximum resident
set size, not allocation size. Efficiency is each version's single-thread
time divided by its parallel time and thread count. Dynamic balancing cuts
16-thread time by 28.1% and 32-thread time by 12.0%; the 8-thread capture is
3.2% slower. It is not a universal speedup or a demonstrated memory reduction.

GNU time average process CPU counts for static/dynamic are 7.64/7.65,
14.57/15.21, and 28.81/28.18 at 8/16/32 threads. These include OpenMP spinning;
they must not be interpreted as useful computational core counts.

## Direct loop evidence

The following occupancy is summed worker loop wall time divided by actual
worker count times the critical loop span, accumulated across calls. It
measures imbalance, not useful CPU utilization or whole-program efficiency.
Timers exclude the following barrier; they include memory stalls and any
preemption inside the work loop.

| Stage | 16-thread static | 16-thread dynamic | 32-thread static | 32-thread dynamic |
| --- | ---: | ---: | ---: | ---: |
| AO-pair contraction | 72.8% | 96.4% | 80.1% | 97.8% |
| Fused AO-H1E action | 70.7% | 95.9% | 82.1% | 97.3% |
| AO-H1E transpose | 71.8% | 95.4% | 83.3% | 97.1% |

At 16 threads, respective phase times change from 9.488/9.638/5.835 seconds
to 6.341/6.097/3.866 seconds. The remaining selected-direction phase grows
from about 0.12 seconds at 8 threads to 2.57 seconds at 32 threads. That is
an identified next target, not addressed by this patch.

## Changes and memory ownership

1. Sort existing unique-integral CSR references in place, then rebuild the
   existing row-reference vector. Integral values are neither copied nor
   duplicated. A standalone 20-action probe (job 280776) reduced average
   primary-edge stride from 25277.7 to 1.99938; sorting cost 0.341 seconds.
   Fused action improvement was about 23%, 15%, and 5% at 8/16/32 threads,
   with relative error below 4.3e-15. This is separate from the dynamic-loop
   comparison above; it must not be counted twice.
2. Partition shared CSR rows into at most workers-squared edge-balanced
   pieces and dynamically dispatch them. Each output row has one owner;
   the existing SIMD contraction and within-row summation order are intact.
   The shared boundary table is about 8 KiB at 32 threads, not a replicated
   integral or operator image.
3. Dispatch scalar, fused, transpose, and block exact AO-H1E integral loops
   with overflow-safe chunks derived from the same piece count. Existing
   AO-sized private accumulators remain; this patch adds no physical-data
   copies. Assignment changes floating-point reduction order, so numerical
   equality is tolerance-based rather than bitwise.
4. Add optional cache-line-sized worker samples under XMVB_PHASE_TIMING=1.
   With timing disabled there are no worker timestamp calls or sample arrays.

The compact non-CSR route and pre-existing block/private workspaces are not
claimed replica-free. Their memory design is outside this measured patch.

## Alternatives checked and rejected

- Jobs 280770 and 280772: compact single-NUMA placement helps somewhat but
  does not eliminate the 16-thread slowdown; therefore affinity alone is
  not the fix.
- Job 280777: PASSIVE wait policy increases 16-thread time from 33.17 to
  38.21 seconds for the original and from 31.44 to 35.89 for sorted traversal.
  No global wait-policy change was retained.
- perf hardware and software events are denied by cluster policy. VTune
  cannot collect under the available ptrace policy. GNU gprofng resolves
  contraction and OpenMP barrier function names, but reports a timer warning;
  its sample percentages are not quantitative time evidence.
- Pending exclusive-node static control job 280799 was cancelled after the
  matched dynamic comparison completed; it supplied no results and would
  not test the latest implementation.

## Verification and reproduction

Job 280804 passes test_active_two_electron_sparse, test_ri_ao_h1e_hvp
(including exact scalar/fused/block actions), and test_orthogonal_direct_ci
at 1, 8, and 32 threads: nine successful executions. Tests cover reference
ordering/orientation, dynamic chunk overflow/small work, and numerical actions.

Release build: -O3 -DNDEBUG -g1 -fno-omit-frame-pointer; Eigen BLAS enabled.
BLAS threads=1, OMP_PLACES=cores, OMP_PROC_BIND=close,
OMP_MAX_ACTIVE_LEVELS=1. No PASSIVE override. CPU sets, executable/input
hashes, output, resource times and phase logs are saved per case.

Remote evidence root:
/home/guqqgroup/taoxia/xmvb-runs/native-20261011/mnf2-storage-280804

Local evidence root: results/studies/locality-20261011/
Parsed final table: balanced.csv. Raw outputs are intentionally not committed.

Use mnf2-storage.sbatch with a source mirror and immutable baseline executable
on Slurm. Its optional `isolated` argument skips wait-policy ablations; it
does not itself request an exclusive node. affinity_probe.py stays inside
the allocated CPU mask. locality_results.py extracts only completed,
converged runs programmatically and computes the comparison CSV.
