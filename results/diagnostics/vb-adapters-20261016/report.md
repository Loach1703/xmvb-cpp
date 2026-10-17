# VB integral adapter cleanup

## Implementation

Removed the exact and RI two-electron builder classes, the VB two-electron
backpropagator and its single-field result (seven files). All consumers now
call the owning common transforms/adjoints. Removed unused injected-component
constructors in the structure, active-gradient and orbital-gradient evaluators,
their unused helper members, two preparation forwarding methods, and duplicate
pair-map/packed-adjoint/RI factor-adjoint implementations, including obsolete
diagnostic copies. No compatibility headers, aliases or fallback were added.

The VB boundary is now `vbscf/integrals/active/space.*`:

- `prepare_timed_active_space_context(input)` performs physical orbital and
  inactive-density preparation, selects exact/RI integrals, and calls common
  transforms directly. It does not thread stateless builder objects through
  multiple evaluators. RI uses the prepared active Eigen block directly;
  its sparse support omits only exact zeros, not a numerical threshold.
- `pullback_active_integrals(input, accepted, packed_gradient)` selects the
  numerical representation and returns an AO-by-active Eigen matrix. It owns
  no duplicate contraction formula or full-AO zero-padded gradient result.
- The common ordinary exact pullback borrows accepted products or consumes
  bounded pair tiles. It builds only the indexing/borrowing metadata needed
  for this operation, not exchange or one-index HVP-only images.
- The ordinary RI pullback uses the existing factor-native response kernel
  instead of the old AO-pair-by-active-pair adjoint and row-major conversions.

The optimizer, coordinates, preconditioner, integral screening, resource
budgets and convergence thresholds were not changed. Contractions now share
the public kernels rather than retaining a separate VB implementation.

## Cluster verification

All compilation and numerical execution used Hanhai22 Slurm compute nodes,
GCC 14.1.0, CMake 3.30.5, strict C++23 and Release. No historical binary was
rerun. Pre-existing local RI worker-budget changes and their additional tests
were excluded and preserved in the working tree.

- Build/regression job **1264365**: exit 0, **41/41 tests passed**, 6.57 s.
- Added ordinary-exact-gradient checks against the accepted HVP cache and
  independently exercised resident and streamed accepted products. Existing
  directional finite differences, gradient/adjoint identities, exact/RI HVPs,
  state-average and molecular optimizer regressions passed unchanged tolerances.
- All ten affected diagnostics compile. Eight were included in the main build;
  the other two compiled in job **1264377**, exit 0. A latent old namespace
  reference in `check_structure_builder_fast_path` was updated to the real
  MCSCF owner, without a compatibility alias.
- **448 production source/header hashes** from the remote build match staged
  Git contents. Standalone common-library and libcint symbol checks passed.

## Small-molecule convergence regressions

Current binary only, eight Slurm-allocated physical cores, OpenMP, one BLAS
thread. Exact integrals, CIAH/Davidson, analytic Hessian-diagonal preconditioner;
projected-gradient threshold 1e-3, energy-change threshold 1e-7 Eh. Original
versioned guesses and basis settings were preserved. No optimizer iteration
limit was introduced by the runner.

| Case | Job | Node | Steps | Orbital HVP | Final energy / Eh | SCF / s | End-to-end / s | Peak RSS / MiB |
|---|---|---|---:|---:|---:|---:|---:|---:|
| 240 | 1264381 | anode23 | 4 | 30 | -343.446302545197 | 8.088943 | 10.963945 | 1528.16 |
| 241 | 1264382 | anode65 | 4 | 44 | -230.720590533578 | 3.103752 | 4.170389 | 904.10 |
| FeCl2 | 1264383 | anode65 | 12 | 220 | -2181.624464845175 | 7.866189 | 8.441548 | 414.46 |

All terminate with `ciah_projected_gradient_tolerance` and exit 0. Initial
energies are -343.075645280883, -230.474471805007 and -2181.617631912268 Eh,
respectively. Saved historical converged energies agree within the existing
1e-7 Eh energy-change scale:

- 240: -343.446302545201 Eh (`orthogonal-tiles-20261015/tile_adjoint_240.csv`, FULL).
- 241: -230.720590533847 Eh (`hanhai22-shared-response-20261012/horizontal-regression-1254764/summary/comparison.csv`).
- FeCl2: -2181.624464866683 Eh (`orthogonal-tiles-20261015/metric_fecl2.csv`).

These are earlier saved revisions, not matched performance baselines; the
timings above are regression observations, not a refactor speedup claim.
The first three launch attempts did not enter SCF: two lacked remote input
decks, and one was rejected by an unnecessary single-LLC pinning constraint.
Inputs were synchronized and the runner now uses Slurm core binding directly.

Executable SHA256:
`27966b5ea621c92392457a56ecd6e25d87167e3f31c56039333b64d963d70b57`.

`evidence.tar.gz` contains build/test logs, source and executable hashes, the
three native inputs/outputs, CIAH traces and GNU time records. It contains no
executables or copied basis library. `tests/vbscf/regression.sbatch` reproduces
the small-system runs using a supplied frozen executable.
