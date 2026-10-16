# C12 transpose support experiments: no retained performance branch

## Decision

Do **not** retain the range/affine experiments as a production optimization.
They pass numerical checks but do not establish an end-to-end speedup over
the saved `0c50888b` parent. All changes to `links.hpp`, `sigma.cpp` and
`sigma_tiles.cpp` were removed. No optimizer, threshold, physical response,
cache budget, tile size, or OpenMP ownership policy was changed.

The retained change strengthens the independent direct-CI regression: shuffled
alpha/beta string enumeration, unequal spin spaces, all exchange parities,
compact sorted output subsets with holes, and chunked nonidentity halos.
`c12_regression.sbatch --tests-only` permits checking this retained test change
without rerunning old end-to-end binaries.

## Workload and controls

Hanhai22 `anode40`, 32 physical cores, affinity `60-91`, OpenMP close/cores,
one BLAS thread. Release C++23 production flags, exact Libcint/cc-pVDZ,
ISCF=9, Davidson. Original input orbitals, thresholds and ITMAX=2000.

Here C12 denotes CAS(12,12) operator fixtures and the actual **LOFLEA HAO**
calculation (132 covalent structures, 924 strings per spin, 150 basis
functions). It is not a separate carbon-C12 molecular geometry. All native
calculations below ran to convergence; they were not one-step probes.
The parent is read from existing job 1264818; no old binary was rerun.

Input SHA256:
`092d545498678b0cf772c5c84b9eab674b694375558cd391a1b739b7bdcc35eb`.
Basis SHA256:
`d53e7e02d55fb99cdf823d50cf13fec375e0cbe1220ff61816aee7ab38e38ec9`.

## Native results

All rows converge in **4 accepted steps, 36 orbital HVPs, 36 structure
responses**. All final energies agree within 3e-12 Eh with the saved parent's
-399.106646267330 Eh. This is numerical equivalence, not a claim of identical
intermediate orbital representatives or bitwise trajectories.

| Variant | Job | Application / s | SCF / s | Peak RSS / MiB |
|---|---:|---:|---:|---:|
| Saved parent | 1264818 | 23.819179 | 21.020280 | 2704.81 |
| Beta range, counting diagnostic | 1264833 | 24.143669 | 21.188099 | 2705.24 |
| Beta range, counters removed | 1264835 | 24.072305 | 21.133800 | 2694.43 |
| Alpha reverse-closure bounds | 1264848 | 24.114454 | 21.257018 | 2725.82 |
| Affine lookup, original ordering | 1264849 | 24.306702 | 21.388568 | 2646.28 |
| Affine lookup, group-count diagnostic | 1264854 | 23.005752 | 20.236521 | 2686.46 |
| Canonicalized complete intervals | 1264855 | 23.862401 | 20.944606 | 2677.07 |

These are distinct successive candidates, not repeated timings of one retained
optimization. Do not select the fastest diagnostic row as a demonstrated
speedup. Native timings and process RSS come from output and GNU `time -v`;
phase runs are separate executions with instrumentation enabled.

## What was strictly removable

1. Beta transpose scatter preserves the alpha spectator. It need not inspect
   rows outside the actual target group's alpha-index bounds. Holes still
   require the exact target map.
2. Alpha transpose scatter changes the alpha index. Its valid bound is the
   reverse single-link closure, not the beta spectator interval. A single
   immutable O(U) bounds array can establish that bound.
3. A consecutive target group has an affine row address, so its explicit
   halo-index map is unnecessary. Only genuinely complete intervals admit
   exact in-place canonicalization to this representation.

The derivation is in `article/channel_write_support.md`. None assumes closure
of the selected VB space, ignores small amplitudes, drops a response, or
replaces H/S by a different operator.

## Why this did not produce the desired gain

The bounded beta counting diagnostic excluded 643255184 of 1353631268 row
candidates (47.52%). Those are *index visits*, not floating-point contractions.
Its atomic diagnostic counters were removed before subsequent production-like
candidates. Removing negative membership lookups does not remove channel
scatter, channel GEMM, Gram products or product-to-structure transforms.

In job 1264854 only 129 of 119196 groups were naturally consecutive. The
canonicalization experiment improves this coverage, but most groups remain
genuinely noncontiguous: only 8127 of 119196 groups (6.82%) become consecutive
in job 1264855. Their maps and their physical contractions cannot be
removed by an interval formula.

All phase runs still execute **129 sweeps, 1171 batches and 1120445436 logical
source-element reads**, with unchanged 1068 product projections. These logical
reads include repeated/structural-zero halo entries; they are not unique
amplitudes or proof of a physical DRAM byte count. The experiments do not
eliminate the dominant source work. This is a source/phase-count conclusion,
not a profiler percentage or a claim that every repeated read is removable.

## Validation and evidence

Each completed candidate passed 9/9 regressions plus selected H/S, structure
response, product-jet finite differences and temporal response checks at
CAS(12,12). Job 1264845 stopped before timing because the memory-layout test
still expected the old link-table byte count; updating that expected metadata
count resolved it in 1264848. It was not a numerical failure.

Raw logs, source/binary hashes, build cache, affinity and numerical output are
archived in `evidence.tar.gz`, excluding binaries. Candidate source snapshots
are archived for reproducibility only and are not compiled or maintained
fallbacks. `comparison.csv` uses the existing native parser;
`phase-summary.json` uses `tests/mcscf/summarize_c12_cost.py`. Nested phase times
overlap and must not be added as independent wall-time shares.

After removing the experimental production branches, a tests-only Slurm job
**1264856** passed **9/9 tests in 3.76 s**. Source hashes of the three restored
production files and both relevant test files match the tested mirror. It
did not rerun an old native SCF calculation. The earlier attempted inline
tests-only submission was rejected by the site's `sbatch --wrap` wrapper
before a job was created; the retained launcher option avoids that wrapper.
