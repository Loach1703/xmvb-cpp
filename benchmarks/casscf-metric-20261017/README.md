# Isolated CASSCF AH metric experiment

Parent: `5c2ad258`. Derivation: `article/casscf_ah_metric.md`.

Outcome: do not switch the production default. See `RESULTS.md`: the exact
AO/MO fixture improves, but the native FeCl2 and F2 tasks regress. To
reconstruct the tested experiment on its parent, apply
`benchmarks/casscf-metric-20261017/metric.patch` with `git apply` from the
repository root. It includes
only the metric integration and its independent tests. There is no retained
alternative/fallback production path. The Hanhai execution mirror and its
`build-final` at this snapshot retain the tested experiment, not the default
local checkout; refresh/rebuild explicitly before a later production run.
Only the generalized AH metric changes from identity to the orthogonal
limit of the existing VB inactive/active span metric: orbital weights
2 for inactive–active and 1 for the other independent rotations; CI weight 1.

Metric images in the common CIAH/coarse workspaces, the shifted diagonal
denominator and the existing keyframe length use this same metric.
Physical gradient/HVP, finite updates, coarse algorithm, recycled seeds,
microstep rules, component cap, coordinate progress norm, outer thresholds,
and convergence criteria remain unchanged. In particular, the .03 rotation
component cap is not multiplied or divided by sqrt(2).

All builds and numerical workloads run in Hanhai25 Slurm allocations.
Reuse `benchmarks/casscf-20261017/build.sbatch SOURCE_DIR BUILD_DIR` for
the qualification build and `sbatch --array=0-2
../casscf-coarse-20261017/run.sbatch SOURCE_DIR` for native F2,
native FeCl2 and the exact FeCl2 AO/MO fixture. Submit from this directory
so artifacts and scheduler logs remain separate. Do not run the baseline
binary again. Reuse its stored records:

- native FeCl2: 291946, 18 macros / 255 HVP, E=-2181.624464862390 Eh,
  solve 3.502806376 s, end-to-end 3.794071006 s, GNU-time wall 3.83 s,
  peak 407.46875 MiB;
- FeCl2 exact fixture: 291944, 16 macros / 209 HVP,
  E=-2181.6244648588327 Eh, solve 12.105322118 s,
  process 12.81 s, peak 1664 MiB;
- native F2: 291953, 3 macros / 14 HVP,
  E=-198.761111541649 Eh, solve .143674790 s,
  process .17 s, peak 56 MiB.

Use the unchanged 32-physical-core C++ environment. Preserve source/binary
and input hashes, node/CPU/affinity, native GNU time reports, convergence
energies and counters. Diagnostic traces are separate from native timings.
PySCF references need not be rerun; their corrected-thread results already
exist in `benchmarks/casscf-coarse-20261017/`.
