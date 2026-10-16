# C12: bounded two-dimensional channel Gram

Theory: [channel Gram tiles](../../../article/channel_gram_tiles.md).
Replace narrow output-row stripes in the exact A^T B cotangent by a two-axis
grid of disjoint output blocks. Each worker reads only the channel columns
needed by its block. No new numerical array, integral replication, larger
workspace, dropped response or changed convergence parameter is introduced.
This reduces redundant logical reads, not the arithmetic asymptotic exponent.

The actual-team-sized grid replaces the old loop; there is no retained
alternative kernel or user tuning parameter. All existing spin multiplicities,
sigma channel reuse and its barrier remain unchanged. Timing-only Gram
barriers are skipped in native runs. Cumulative scope times are nested.

## Qualification

Hanhai22 Slurm **1264607**, `anode40`, 32 physical cores, affinity `60-91`,
GCC14.1, C++23 Release with the existing optimization flags, BLAS threads=1.
Completed exit0. **9/9 CTests passed** (3.26s), including equal/unequal spins,
exchange parities, multistate responses, nonclosed selected structures and
narrow evicting workspaces. Independent CAS12 checks:

| Check | Error norm |
|---|---:|
| Selected H action | 1.97942125e-14 |
| Selected M action | 1.90486032e-15 |
| Structure-response image | 3.36632466e-16 |
| Structure-response integral adjoint | 1.47706187e-15 |
| Orbital HVP vs gradient finite difference | 5.84911563e-11 |
| Lagrangian directional derivative | 5.21188648e-12 |

The optional full-carrier product-jet replay was disabled; its printed zero
placeholder errors are not used as evidence. Independent finite differences
and separate H/M/response references execute and pass.

## Native end-to-end convergence

Same byte-identical original LOFLEA CAS(12,12), COV132, HAO, exact cc-pVDZ,
ISCF=9 and Davidson. Same tolerances and complete convergence, not a truncated
SCF. Only new versions were run: baseline data below were already saved.

| Metric | Start: 6a7b74da / 1264598 | Direction reuse: 977da393 / 1264605 | New Gram / 1264607 |
|---|---:|---:|---:|
| Accepted steps | 4 | 4 | 4 |
| Orbital HVP / physical CI responses | 36 / 36 | 36 / 36 | 36 / 36 |
| Final energy / Eh | -399.106646267330 | -399.106646267331 | -399.106646267331 |
| SCF wall / s | 32.038107 | 31.425332 | 29.927104 |
| Application end-to-end / s | 34.899360 | 34.213118 | 33.073310 |
| Process peak RSS / MiB | 2697.812 | 2663.473 | 2723.004 |
| Mean accounted cores | 27.67 | 27.55 | 27.14 |

The Gram change reduces native application wall by 3.33%, SCF wall by 4.77%.
Both changes together reduce native application wall by **5.23%**, SCF wall
by **6.59%**. Final energies differ by about 1e-12 Eh. The per-step HVP counts
are unchanged at 6,17,10,3. Peak RSS is 25.19MiB (0.93%) above the starting
run, so there is no claimed RSS reduction; no new explicit numerical buffer
is allocated by the Gram change. BLAS views/packing and process allocator
peaks are distinct from the explicit shared workspace bound.

## Separate phase-enabled attribution

| Inclusive scope / s | Saved direction reuse | New Gram |
|---|---:|---:|
| 36 complete coupled actions | 15.721623 | 13.377332 |
| Physical orbital HVP | 9.634859 | 8.465292 |
| Structure-side scope | 6.065116 | 4.887558 |
| 88 sigma/integral-adjoint sweeps | 13.920101 | 10.910750 |
| Tile compute, including Gram | 9.573236 | 6.364395 |
| Tile read | 3.908233 | 3.695451 |
| Tile write | 2.485603 | 2.410011 |

Tile compute falls 33.5%, sigma/adjoint scope 21.6%. Ordinary H/S sigma and
physical-response sweep counts are unchanged: 165 total sweeps, 3109 batches
and 1,964,120,314 logical source elements. The new measured Gram scope is
**2.954028s**, nested inside tile compute, with 114176 output-block tasks.
Older data do not have a Gram scope; it is missing, not zero. No estimate of
measured DRAM traffic or sampling-profiler percentage is asserted.

Remaining tile read costs about 3.70s and merits the next review: retained
directional input fits on the same support within the accepted two-component
budget, but its reader still schedules separate copies per source panel.
Check whether these can be borrowed and copied in one batch without new
numerical storage. Keep partial/evicted lifetime constraints exact.

`comparison.csv`, `phase-summary.json` and `evidence.tar.gz` retain parsed and
raw results, input/binary/source hashes, environment, build and CTest logs.
No old source/executable was rerun and no unrelated local RI edits were synced.
