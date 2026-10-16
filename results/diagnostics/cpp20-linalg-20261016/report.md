# C++20, CIAH naming, and shared numerical modules

Parent: `385566af`. This refactor changes interfaces and ownership, not the
orbital coordinates, derivative formulas, optimization decisions or tolerances.

## Interfaces

- Require C++20, disable language extensions, and remove `-fpermissive`.
  The tested GCC14.1 build uses `-std=c++20`, `-O3 -DNDEBUG`, `-g1` and
  `-fno-omit-frame-pointer`; there is no numerical flag relaxation.
- ISCF=9 is **CIAH**, with `--optimizer-backend ciah`, `--ciah-trace`,
  `CiahIterationRecord`, `CiahTraceWriter` and `XMVB_CIAH_TRACE`. Original
  CLI spellings are rejected, not retained as aliases. Numeric ISCF values
  remain 5=L-BFGS, 7=TNHVP, 8=P-LBFGS, 9=CIAH.
- ISCF=8 is **P-LBFGS**, preconditioned L-BFGS, CLI `plbfgs`. Its local
  orbital-curvature initial inverse and transported secant updates are
  unchanged. `BlockInverseBfgs` is a different, genuinely block-secant
  operation and retains its descriptive name.
- `QuadraticProblem` is shared by CIAH and the distinct norm-extended solver.
  CIAH lives in `mcscf/optimization/ciah/solver.*`; the actual NEO
  trust-region solver remains in `mcscf/optimization/neo/`.
- `xmvb_linalg` owns generalized dense/Davidson eigensolvers, eigen response,
  metric solves, MINRES, orthonormal bases and projection symmetry bounds.
  `linalg::VectorAction` is the shared numerical callback. `linalg/` and
  `parallel/` cannot depend on an electronic-model layer; CMake enforces this.
- OpenMP execution helpers moved unchanged to `parallel/openmp.hpp`.
  No new MPI path or worker-local large arrays were introduced.
- Dense eigensolver input is a read-only `std::span<const double>`; vectors,
  arrays and other contiguous storage can be borrowed without prescribing
  ownership. LAPACK's required writable H/M copies are unchanged. Eigen
  remains the dense physical matrix/vector representation.

The structure and naming contracts are in `docs/code_organization.md`.
Manuscript drafts and saved benchmark labels are not rewritten retroactively.
For the language/toolchain choice, see the official
[GCC support table](https://gcc.gnu.org/projects/cxx-status.html) and
[CMake standard selection](https://cmake.org/cmake/help/latest/prop_tgt/CXX_STANDARD.html).

## Build and regressions

Hanhai22 final job **1264177**, anode42: **33/33 tests passed**. Main program
and both changed diagnostic tools build without permissive C++ compilation.
New tests cover span input/extent checks and shared-buffer execution inside
an existing OpenMP team. Existing regressions cover the ISCF mapping, accepted
and rejected CLI names, public solvers, Davidson/eigen response, CASCI,
transport, structure tiles, F2 HAO/OEO CIAH, matched-seed PySCF-reference F2,
state-average and text/trace/orbital output.

`test_minres`, `test_davidson` and `test_eigen_response` link only the numerical
library. Their saved symbol tables have zero VB/MCSCF/CASSCF function
definitions. CIAH/NEO/CASCI tests also contain no VB definitions. Source
hashes match497 local files. The only two differences are the pre-existing
local RI changes; just their necessary OpenMP include-path edits are staged.
The original RI implementations, not unrelated local modifications, were
used in this cluster qualification.

Intermediate build issues were new GNU-make goals requiring explicit CMake
configuration and a double-qualified test helper after namespace migration.
Both were corrected. The execution mirror also lacked the versioned F2_SA2
fixture; transferring the repository's existing input fixed that test setup.
No numerical threshold or expected energy was changed to pass a test.

## Native end-to-end COV convergence

LOFLEA CAS(12,12), HAO,132 COV structures, original input orbitals, exact
libcint cc-pVDZ, CIAH, Davidson and analytic Hessian-diagonal. Eight pinned
OpenMP cores in one LLC, BLAS one thread. Gradient tolerance1e-3 and
energy-change tolerance1e-7 Eh. Saved versions were not executed again.

| Version / job | Node | Steps | HVP / response | Final energy (Eh) | Process (s) | RSS (GiB) |
|---|---|---:|---:|---:|---:|---:|
| Before this refactor /1264100 | anode42 | 4 | 36 /36 | -399.106646267329 | 51.26 | 2.508 |
| C++20/shared modules /1264182 | anode42 | 4 | 36 /36 | -399.106646267332 | 48.17 | 2.495 |

Final energy agrees within3e-12 Eh. Per-macro microsteps (10,31,17,3), HVPs
(6,17,10,3), keyframes (3,8,3,1) and rejected trials (1,1,0,0) are unchanged.
This one native regression shows no wall-time or memory regression; it is not
evidence for a general speedup from the language-standard change. GNU-time
peak RSS is reported above; Slurm polling is retained separately. SCF wall
is44.217137 s, application wall47.976280 s, process wall48.17 s.

Input SHA-256:
`5d3086bf489db0d039c33197dc3ca61f3e5dc85bc4e1a1aaada292383ef0a858`.
Executable SHA-256:
`6b9b5c966a16c5da1cd0aa73a14d62007f4bb0ea577cdaeaf6a3ee9b416eab4f`.
The evidence archive contains raw build/test logs, compiler flags, numerical
library symbol tables, source/binary hashes, native input/output/trace,
affinity and GNU-time data. `native.csv` uses the existing automated parser;
the trace filename `neo.tsv` is an arbitrary historical benchmark filename,
not an accepted old CLI spelling.
