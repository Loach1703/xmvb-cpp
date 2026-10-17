# Exchange-reduced integral adjoint and complete response

## Implementation

Parent revision: `4de4db09`. No optimizer, coordinate, tolerance, or integral
screening changes. The new work is in `mcscf/ci/sigma.cpp` and its existing
structure-action consumers.

- Extend the exact spin-exchange quotient to the integral adjoint and fused
  response. The structure expansion certifies the common parity once.
- Pack only valid triangular rows into shared tiles, using arithmetic row
  offsets. No full carrier index array, zero-filled boundary work, or new
  per-thread matrix is required.
- Restore the **complete** response image with the exact sector projector.
  The representative used for a final structure projection is not itself
  the physical response; this distinction is tested explicitly.
- Fuse accepted-point sigma and the integral adjoint. Identical operands
  share one link scatter instead of a standalone sigma and two adjoint
  scatters. Unequal spin spaces also benefit from this fusion.
- Retain the general contraction for mixed sectors and unequal spin lists,
  where exchange reduction is not a valid identity.

Theory: [spin_exchange_contraction.md](../../../article/spin_exchange_contraction.md).

This reduces a constant factor, not a scaling exponent. For CAS(12,12),
924 strings per spin give 853776 product rows; the even/odd representative
counts are 427350/426426. The combined shared tile bound remains 64 MiB.

## Verification

Final packed-row candidate: Hanhai22 job **1260303**, anode26, 8 bound cores,
GCC14.1.0, Release `-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`, BLAS1.
Build and all calculations ran under Slurm. Sources were placed in a new
mirror so jobs1260267/1260268 still refer to the preceding implementation.

All eight CTests passed: product_frontier, ci_transport,
representation_overlap, restricted_metric, selected_ci_transport,
orthogonal_direct_ci, structure_tiles, coupled_structure.

The new tests compare both integral adjoints and the entire sigma against
unreduced current-code contractions, including a separate ordinary sigma.
They cover both parities, nonidentical left/right inputs, aliasing accepted
inputs, multiple RHS, strided views, and partial tiles. Existing tests cover
mixed selected spaces, unequal spins, independent derivatives, and coupled
structure responses.

Maximum absolute difference across the benchmark's sigma/integral-adjoint
comparisons is 2.23e-16; for CAS(12,12) it is 1.91e-17. Checksums include
the complete image, not just a scalar selected-space projection.

## Native kernel timings

The table compares the **same new binary** with and without a certified
exchange sector. The unreduced option is the general mathematical domain,
not a rerun of a historical version. Both accepted columns already use
scatter reuse; the table isolates quotient savings, not the separate gain
from accepted-point fusion. Each entry is one warmed call on anode26.

| Orbitals | RHS | Parity | Operation | Full domain, ms | Quotient, ms |
|---:|---:|:---:|---|---:|---:|
| 10 | 2 | + | response sigma + adjoint | 14.588 | 8.696 |
| 10 | 2 | + | accepted sigma + adjoint | 13.673 | 8.006 |
| 12 | 1 | + | response sigma + adjoint | 238.412 | 166.562 |
| 12 | 1 | + | accepted sigma + adjoint | 216.624 | 133.507 |
| 12 | 1 | - | response sigma + adjoint | 275.897 | 139.129 |
| 12 | 1 | - | accepted sigma + adjoint | 218.942 | 111.940 |

The full raw CSV also contains the smaller cases and the odd 10-orbital
measurements; no slower sample is removed. These are short kernel timings,
not molecular SCF runtimes. Whole-process peak RSS is **179.81 MiB** for the
combined validation/benchmark process containing both variants and reference
arrays; it cannot establish an isolated per-variant RSS change.

Preliminary rectangular-boundary candidate1260298 also passed all eight
tests; it was superseded by packed boundary rows, not retained as another
production algorithm.

## End-to-end status

Job **1260305** is queued for the same anode15/32-core setup as saved
240 job1260261. It will run COV, ION(0:1), ION(0:2), ION(0:3), and FULL
to convergence and record energy, steps, HVP, native time, and RSS.
Historical binaries are not rerun.

No new end-to-end speedup is claimed yet. Previous phase logs show that
CI integral-adjoint/response time is only tens of milliseconds for 240;
this change is primarily useful to larger carrier contractions.
