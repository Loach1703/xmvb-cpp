# FeCl2 PySCF thread diagnosis, 2026-10-17

The existing full PySCF FeCl2 jobs allocated 32 cores, but GNU time reported
about one effective core: coupled Newton job 291883 took 1194.60 s wall,
411.52 s user + 782.95 s system, 823692 KiB peak RSS, and converged to
-2181.624464854347 Eh (24 macro callbacks, cumulative `tot_hop=326`).
Standard AH job 291887 took 1665.37 s wall, 575.38 s user + 1089.75 s
system, 851976 KiB peak RSS, and converged to -2181.6244648112297 Eh
(final macro index 35; 121 callback entries). Those full jobs used
`OMP_NUM_THREADS=32`, `OMP_PROC_BIND=close`, `OMP_PLACES=cores`, and serial
MKL/OpenBLAS environment settings. Their separate `affinity.txt` files show
the *probe process* had a 32-CPU affinity, not the later Python process.

## Bounded diagnostic

Two 32-core/96-GiB Slurm allocations on `anode029` reused the exact frozen
FeCl2 AO/MO fixture and original PySCF 2.14.0 script. The diagnostic wraps
only one standard AH macro iteration, with a fresh Python process per variant.
It intentionally stops before convergence. All variants report the same
one-macro energy to within 5e-12 Eh: -2181.61762853247 Eh; all have four
callback entries. The PySCF wheel's OpenBLAS reports
`threading_layer=disabled` and stays at one thread. NumPy/SciPy load Anaconda
MKL; both Intel and GNU OpenMP runtimes are present.

| Job / variant | OMP | MKL | Binding | AH kernel wall (s) | Process CPU (s) | GNU-time whole-process wall (s) | Peak RSS (KiB) |
| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| 291921 / BLAS1 | 1 | 1 | close + cores | 8.517 | 8.458 | 10.28 | 811424 |
| 291921 / BLAS32 | 1 | 32 | close + cores | 49.992 | 49.989 | 51.56 | 817628 |
| 291926 / original-style control | 32 | 1 | close + cores | 84.386 | 84.329 | 86.36 | 839844 |
| 291926 / unbound OpenMP | 32 | 1 | false, no PLACES | 1.337 | 36.333 | 2.03 | 791164 |
| 291926 / unbound MKL | 1 | 32 | false, no PLACES | 10.263 | 252.876 | 10.93 | 780044 |

The decisive matched comparison is the two OMP=32/MKL=1 variants in job
291926. Before NumPy/PySCF import, the process had a 32-CPU affinity in both.
With `OMP_PROC_BIND=close` and `OMP_PLACES=cores`, its main thread then had
one CPU and repeated `/proc/self/task` samples found all 33 threads permitted
only CPU 120. GNU time reported 99% CPU. With `OMP_PROC_BIND=false` and
`OMP_PLACES` unset, the main thread and all 33 sampled threads retained masks
containing the 32 allocated CPUs; GNU time reported 2025% CPU. Both variants
used the same `srun --cpu-bind=cores`, PySCF `lib.num_threads(32)`, fixture,
source, MKL=1, and wheel. The binding-setting pair is therefore the measured
cause of single-core effective execution in this representative phase. The
experiment does not isolate which of the two binding variables is responsible,
or whether the coexistence of Intel/GNU OpenMP runtimes contributes.

The job 291921 BLAS32 variant changed MKL to 32 but left the PySCF bundled
OpenBLAS at one thread; it still ran near one-core CPU usage under the binding
settings. In job 291926, unbound MKL32 used about 23.6 cores on average but
was slower than unbound OMP32/MKL1 for this phase. Thus setting all BLAS to
32 threads is not a supported remedy.

For this PySCF runtime and fixture, use
`OMP_NUM_THREADS=32 OMP_PROC_BIND=false OMP_MAX_ACTIVE_LEVELS=1`
with `OMP_PLACES` unset and
`MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1`.
Keep `lib.num_threads(32)` and the 32-core Slurm request. The 1.337 s result
is one AH macro only; no full solver convergence or end-to-end speedup was
measured under the corrected environment.

Raw records are in `results/291921` and `results/291926` (wheel extraction
excluded from the local copy), and on Hanhai25 at
`/home/guqqgroup/taoxia/work/xmvb-ciah-subspace-20261017/benchmarks/casscf-threads-20261017/pyscf-blas-291921`
and the adjacent `pyscf-affinity-291926`. Each has script, fixture, metadata,
and wheel SHA-256 values; per-variant logs, GNU time, threadpool/library
metadata, and measured per-thread affinity samples. The Slurm mirror's
`parent-revision.txt` marker was `146dfc87`, also recorded in the earlier
full PySCF job 291883 artifact. The diagnostic jobs themselves recorded
script, fixture, metadata, and wheel hashes, but did not snapshot a C++ source
revision. This mirror marker is separate from the current C++ implementation
checkout at `2797d87b` (with additional uncommitted changes), and from the
earlier qualification report's `65c1320b` parent. The fixture metadata has no
source revision field. PySCF wheel SHA-256 was
`37b0bccc55450311a55318cd643e851353331ddeab4fc0c0065e83c905e41502fa57a`.
Slurm jobs 291921 and 291926 both completed successfully; no profiler was
needed because native timing and direct affinity/threadpool observations
resolved the question.
