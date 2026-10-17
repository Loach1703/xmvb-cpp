# Common integral kernel extraction

## Changes

- One `xmvb_integrals` library owns AO exact/RI actions, AO-pair graph and
  backend contracts, active one-/two-electron transforms, integral directions,
  adjoints and fixed-adjoint HVPs. Public namespace: `xmvb::integrals`.
- Exact sparse support is borrowed through `std::span`; RI coefficients use
  Eigen matrices. The kernels do not require orthonormal orbitals or VB input
  objects. No change to integral screening, tiling, memory budgets, precision,
  preconditioning or optimizer decisions was made.
- Exact/RI forward transformation implementations were removed from VB builders;
  those adapters supply VB orbital-preparation data to the common functions.
- `xmvb_libcint` now links the integral library, not `xmvb_vbscf`. Backend APIs
  use `xmvb::libcint`. Generic storage conversion and worker timing belong to
  `linalg/storage.hpp` and `parallel/workers.hpp` respectively.
- Old moved paths were deleted, without forwarding headers or namespace aliases.
  Two unused allocating overloads and a test-only pair-map forwarding function
  were deleted. The test calls the owning pair-map implementation directly.
- The AO-to-CASCI preparation entry uses the shared AO response and active
  transforms. This is not a completed CASSCF orbital optimizer.

## Verification

All compilation and numerical execution used Hanhai22 Slurm compute nodes,
GCC 14.1.0, CMake 3.30.5, strict C++23, Release, OpenMP and BLAS/LAPACKE.
No old binary was rerun.

- Final build/regression job: **1264330**, node **anode56**, exit 0.
- **41/41 tests passed** in 6.36 s.
- Shared tests cover exact/RI forward transforms, integral directional finite
  differences, adjoint duality, fixed-adjoint HVP finite differences and
  exact/RI HVP agreement. CASCI preparation is compared with an independently
  transformed nontrivial orthonormal MO frame, including inactive/nuclear energy.
- Existing VB tests cover HAO/OEO diagonal preconditioning, selected-structure
  actions/transport, sparse exact transforms, RI responses, state averages,
  Davidson and CIAH/TNHVP/P-LBFGS molecular regressions.
- `nm -C` checks passed for standalone integral, libcint, CASCI and shared solver
  tests: no defined `xmvb::vb` implementation symbols.
- All **455 production C++ source/header hashes** in the remote build manifest
  match the staged Git contents. Pre-existing local RI worker-budget changes
  and their test additions were excluded; they remain in the working tree.

An expanded old test initially assumed every resident batch must return a pair
image. The existing accepted-exchange path instead directly returns a fixed
adjoint. The test now explicitly checks both contracts against the scalar HVP,
including a forced resident pair-image case. Production algorithms and test
tolerances were not changed to make that test pass.

## End-to-end regression

Job **1264335**, anode56, 8 pinned physical cores, one BLAS thread. LOFLEA
CAS(12,12), HAO, 132 covalent VB structures, exact cc-pVDZ integrals,
CIAH/Davidson, projected gradient 1e-3 and energy change 1e-7 Eh.

| Quantity | Current extraction | Saved job 1264182 |
|---|---:|---:|
| Accepted steps | 4 | 4 |
| Orbital HVP actions | 36 | 36 |
| Structure response actions | 36 | 36 |
| Final energy / Eh | -399.106646267330 | -399.106646267332 |
| SCF wall / s | 50.305233 | 44.217137 |
| Application wall / s | 54.394423 | 47.976280 |
| Process wall / s | 54.63 | 48.17 |
| Peak RSS / KiB | 2638048 | 2615788 |
| Node | anode56 | anode42 |

The saved run used the preceding C++20/linalg revision, not the current C++23
binary. Different nodes and runs prevent attributing the timing difference to
this refactor; this is a correctness regression, not a speedup qualification.
The final energy difference is approximately 2e-12 Eh. Per-macro microiteration,
keyframe, HVP and rejection counts are unchanged:

| Macro | Microiterations | Keyframes | HVP | Rejections |
|---|---:|---:|---:|---:|
| 1 | 10 | 3 | 6 | 1 |
| 2 | 31 | 8 | 17 | 1 |
| 3 | 17 | 3 | 10 | 0 |
| 4 | 3 | 1 | 3 | 0 |

Executable SHA256:
`1b8840c56e8dad36c7b14a8e4015e722adbcd1702a60214cfd954a5ee29df108`.

`evidence.tar.gz` contains build/test logs, source/binary hashes, and the current
COV input, output, CIAH trace, resource record and affinity information. It does
not contain executables. The saved reference is recorded in
`results/diagnostics/cpp20-linalg-20261016/native.csv`.
