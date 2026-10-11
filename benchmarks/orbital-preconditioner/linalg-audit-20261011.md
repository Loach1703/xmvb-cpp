# Eigen and BLAS/LAPACK audit

## Decision

Keep Eigen containers, views, and expression interfaces. Choose explicit
BLAS/LAPACK routines only where measured algorithm, workspace, or data-layout
differences justify doing so. Do not replace every Eigen expression.
This round adds opt-in scalar shape counters, not numerical-kernel changes.

Hanhai25 builds enable `EIGEN_USE_BLAS` and `EIGEN_USE_LAPACKE` and link
OpenBLAS. The installed Eigen source maps dynamic double LLT to `dpotrf`,
SelfAdjointEigenSolver to `dsyev`, JacobiSVD to `gesvd`, and BDCSVD to `gesdd`.
Executable imports confirm BLAS GEMM/GEMV and LAPACKE factorization symbols;
imports alone do not establish execution counts. These mappings are version,
scalar type, and expression dependent; small or special cases can bypass them.

Reference: [Eigen backend documentation](https://libeigen.gitlab.io/eigen/docs-master/TopicUsingBlasLapack.html).

## Measured MnF2 inventory

Job `281282` ran exact Libcint, NEO/Davidson, on `anode020` with 8, 16,
and 32 cores. `XMVB_MATRIX` records operation-site counts and shape extrema,
not intercepted LAPACK calls. All three runs have the same shape counts.
Row and column extrema of rectangular operations are independent ranges.

| Operation | Sites | Observed dimensions | Backend in this build | Assessment |
| --- | ---: | --- | --- | --- |
| NEO relaxed spectrum | 694 | square, 1 through 231 | explicit `dsyevd` | Highest-priority dense spectral site; previous matched run established a gain over `dsyev` |
| NEO metric factor | 669 | square, 1 through 231 | Eigen LLT / `dpotrf` | Already LAPACK-backed; bypassing Eigen is not by itself an improvement |
| NEO structure pseudoinverse spectrum | 39 | square, 0 through 18 | Eigen / `dsyev` for nontrivial sizes | Tiny; do not prioritize replacement |
| Strict-locality constraint SVD | 775 | rows 72-88, columns 17-18 | Eigen JacobiSVD / `gesvd` | Audit rank stability, not just timing |
| Local gauge complement SVD | 775 | rows 4-10, columns 36-52 | Eigen JacobiSVD / `gesvd` | Full right vectors are required for the nullspace; thin vectors cannot replace them |

Locations: `optimization/neo/response_solver.cpp` and
`orbitals/charts/chart.cpp`. Only these five operation sites have new shape
counters. Counts for other factorization sites remain unmeasured; source
occurrences must not be presented as dynamic invocation counts.

At maximum size, one 231-by-231 double matrix has about 0.407 MiB of data.
That is not the eigensolver's total workspace or process peak RSS. The shape
counter stores no matrices and introduces no per-thread numerical buffers.

The previous matched benchmark (`281273`, source `dba4c99`) found the 32-core
projection-factor stage reduced from 1.436 to 0.652 s. That includes reuse
and driver changes; it is not a pure Eigen-versus-LAPACK comparison.
See `neo-projection-reuse-20261011.md` for the developmental comparisons.

## Other production paths

| Path | Present representation / backend | What to investigate |
| --- | --- | --- |
| Direct-CI sigma | Eigen `noalias()` products, configured BLAS; bounded tiled scatter/reduction | Shapes, small-GEMM packing, team granularity and reuse; do not assume rewriting `dgemm` call syntax speeds it up |
| Direct-CI integral adjoint | Shared rectangular tiles and disjoint output blocks; BLAS-compatible products | Actual strides, packing and tile utilization; preserve the shared-memory budget |
| Physical metric inverse | Certified iterative CG and cached block preconditioner | Repeated operator/projection work, not wholesale conversion to a dense LAPACK solve |
| Exact AO-pair action | Sparse graph rows and SIMD kernel | Not a dense GEMM problem; densifying the graph increases memory and discards sparsity |
| AO one-electron actions | Specialized fused/transpose sparse-integral loops | SIMD and stream reuse, not a generic LAPACK replacement |
| Block curvature / shifted preconditioners | Eigen LDLT / LLT and small spectra | First measure dimensions and factor counts; keep positivity and rank semantics |
| Spin-pair overlap rank resolution | FullPivLU and BDCSVD | Preserve full-pivot rank/determinant/nullity conventions; LAPACK partial-pivot LU is not an interchangeable substitute |

Core Davidson already explicitly selects projected low roots with `dsyevr`
in `core/eigensolver.cpp`. The RI metric uses explicit `dsyevd` in
`libcint/ri_provider.cpp`. There is no need to invent another eigensolver
framework for these existing operations.

The spin-pair overlap row is a source inventory entry, not evidence that it
executes on MnF2's orthogonal direct-CI path. Its dynamic counts were not
collected in this run.

Eigen LDLT does not have the same pivot semantics as a general LAPACK
symmetric-indefinite factorization. A positive-definite path can be evaluated
with Cholesky only after proving and retaining its required checks. Likewise,
the documented LAPACKE substitution for JacobiSVD changes the numerical
algorithm; an accuracy audit must include near-rank-deficient matrices,
nullspace residuals, dimensions, and downstream gradients, not just energies.
No SVD or rank threshold was changed in this audit.

## Prioritized experiments

1. Measure allocation/workspace cost at the 694 NEO spectrum sites. Test
   caller-owned LAPACK `_work` buffers on observed subspace sizes, retaining
   the same driver and full residual certification. Existing phase data do
   not yet establish that allocation is a dominant cost.
2. Test the two strict-locality SVD sites for rank and nullspace stability
   around the existing thresholds. If a more robust algorithm is required,
   make an explicit numerical choice, not a global speed-only backend switch.
3. Collect GEMM dimensions/strides on a larger carrier before changing
   direct-CI BLAS calls. MnF2 alone cannot justify a general tile planner.
4. Only then consolidate proven repeated operations into minimal core helpers
   taking Eigen storage/views and explicit workspace ownership. Keep full
   spectra and selected roots as distinct operations. Avoid a backend factory,
   hidden retries, molecule-specific cutoffs, or additional full-array copies.

## Audit run and reproduction

All runs converged in 14 steps with 669 orbital HVPs, at unchanged energy and
gradient tolerances of 1e-7 hartree and 1e-3.

| Cores | Instrumented wall, s | Energy, hartree | GNU-time peak RSS, MiB |
| --- | ---: | ---: | ---: |
| 8 | 25.877223 | -1348.893353119307 | 453.2 |
| 16 | 16.303401 | -1348.893353119321 | 449.4 |
| 32 | 12.607264 | -1348.893353119316 | 489.3 |

These are audit runs, not a new optimization speedup claim. The node was not
exclusive; compact affinity spans one NUMA domain at 8/16 cores and two at
32 cores. BLAS uses one thread to avoid nested parallelism. The sparse exact
two-electron and RI AO-H1e HVP tests passed in the compute-node batch.

Run `mnf2-balance.sbatch SOURCE_DIR`, then parse logs with `linalg_results.py`,
`phase_results.py`, and `locality_results.py`. Raw local evidence is under
`results/studies/locality-20261011/mnf2-balance-281282/`; shape CSV is
`linalg-281282.csv`. Remote evidence is under
`/home/guqqgroup/taoxia/xmvb-runs/native-20261011/`.
