# C12: single-write halo readers and fused RHS packing

## Scope

Hanhai22 `anode40`, 32 physical cores, CPU affinity `60-91`, OpenMP
`close/cores`, BLAS threads=1. GCC14.1/C++23 Release production flags.
The saved parent is `0c50888b`, job 1264818. It was **not rerun**.

The operator fixtures are CAS(12,12), COV132. The actual converged molecular
calculation is **LOFLEA HAO**, exact Libcint/cc-pVDZ, 150 basis functions,
132 covalent structures, 924 strings per spin, ISCF=9/Davidson. This is not
a separate carbon-C12 molecular input. Original orbitals, convergence
thresholds and ITMAX=2000 remain unchanged. Numerical runs and builds use
Slurm compute nodes only. Unrelated local RI changes were not transferred.

Input SHA256:
`092d545498678b0cf772c5c84b9eab674b694375558cd391a1b739b7bdcc35eb`.
Basis SHA256:
`d53e7e02d55fb99cdf823d50cf13fec375e0cbe1220ff61816aee7ab38e38ec9`.

## Exact dataflow argument

The sigma halo reader must produce one accepted-state RHS lane and literal
zeros in unused lanes. All implementations of the local `y` reader fully
overwrite their requested output: full/resident columns are copied, compact
columns clear their structural holes before scattering, and unreachable
columns are explicitly zero. Therefore clearing the **selected** lane before
calling `y` is a redundant second write. A one-lane request needs no wrapper
clearing at all; a two-lane request clears only its unused lane. The CI-only
reader uses the same borrowed adapter, without copying numerical owners.

The generated primal/orbital/CI RHS panels used to be concatenated into a
second dense matrix, then gathered into their exact compact support. Selection
acts independently on each RHS, so the common gather can instead read each
original panel directly. This removes the concatenated dense owner, its
initialization and its copy. The one-matrix full-support move remains a
zero-copy path. Both overloads use one numerical gather implementation.

Two further lifetimes are unnecessary: cached calls never need to reconstruct
the primal panel, and accepted-only calls never consume a tangent panel. They
are no longer allocated. When an anchored orbital tangent consists solely of
the one-body generator image (with the CI component either absent or separate),
the generator fully overwrites the tangent directly. No zero tangent plus
separate generator owner is needed. A genuinely additive orbital+CI tangent
still retains its required summation.

No physical overlap, H/S, integral adjoint, orbital response, CI response,
transport term or exchange component is removed. No amplitude or residual
screening, selected-space closure assumption, new U-squared cache, thread-local
carrier, larger tile, precision change or optimizer-policy change is involved.
The existing shared workspace budgets and task ownership are unchanged.
This is a memory-traffic/dataflow improvement, not a new asymptotic scaling
claim or proof of improved parallel scaling.

## Stages and evidence

- **1264867:** selected-lane zeroing removed, shared interleaved reader.
- **1264868:** no initialization of a fully overwritten concatenation buffer;
  omit unused primal/absent tangent panels. This intermediate candidate still
  contains the concatenation buffer.
- **1264871:** final direct multi-RHS gather, generator writes its actual
  destination, no concatenated dense RHS owner.

Native timings and process peak RSS are parsed with the existing native
parser; separate phase-enabled executions locate changed scopes. Nested
phase scopes overlap and must not be summed as independent wall-time shares.
Logical source-element counts are not physical DRAM traffic. No sampling
profiler percentage or hardware-bandwidth claim is made.

All results and final validation are recorded in the generated CSV/JSON and
the qualification section below. Raw evidence includes source/binary/input
hashes, compiler cache, affinity, unit logs, CAS12 operator checks and both
native/phase SCF runs; binaries are excluded from the archive.

## Qualification

Job **1264871** completed exit0. **9/9 regressions passed in 3.64 s**. Direct
multi-panel packing also has a new exact-equality test for compact rows,
empty columns, tiny entries (1e-180), multiple RHSs and full support.
The existing shuffled enumeration, compact output, unequal-spin, exchange,
small-workspace and multi-state regressions remain active.

The final source hashes of `product_jet.cpp`, `column_tiles.hpp` and the
modified unit test match the compiled mirror manifest. Independent CAS12
checks retain the following errors (norms, not energy tolerances):

- Physical structure-response image: 3.4002682088792616e-16.
- Structure-response integral adjoint: 1.2707391900175511e-15.
- Orbital HVP versus gradient finite difference: 2.5933568576767072e-11.
- Lagrangian finite difference: 6.9466100081384818e-12.

The optional complete-carrier jet replay is not enabled in this launcher;
its placeholder zeros are not additional parity measurements. The independent
finite-difference and structure-response reference checks are actually run.

## Native converged comparison

| Stage | Job | Steps / HVP | Final energy / Eh | SCF / s | Application / s | Peak RSS / MiB |
|---|---:|---:|---:|---:|---:|---:|
| Saved parent | 1264818 | 4 / 36 | -399.106646267330 | 21.020280 | 23.819179 | 2704.81 |
| Single-write reader | 1264867 | 4 / 36 | -399.106646267330 | 20.209089 | 22.874827 | 2702.17 |
| Omit unused panels | 1264868 | 4 / 36 | -399.106646267329 | 20.379011 | 23.310588 | 2673.28 |
| Direct lane gather | 1264871 | 4 / 36 | -399.106646267329 | 19.845801 | 22.626016 | 2687.16 |

Every run uses 36 physical structure responses, per-step HVP counts 6,17,10,3,
and the original `ciah_projected_gradient_tolerance` termination. Final native
application time decreases **5.01%**, SCF time **5.59%** versus the saved parent.
Peak RSS is essentially unchanged; the smaller staging lifetime is not a
claim of a substantial measured whole-process RSS reduction. These are single
matched measurements, not a statistical confidence interval. No old version
was rerun. Intermediate candidates are shown rather than selecting only the
fastest row.

## Separate phase-enabled comparison

| Nested scope / s | Saved parent | Single-write reader | Final |
|---|---:|---:|---:|
| Application | 24.085821 | 22.701306 | 22.605579 |
| Complete coupled actions | 7.599649 | 6.589541 | 6.522412 |
| Physical orbital HVP | 7.533450 | 6.526179 | 6.458665 |
| Integral-direction sweep | 6.443392 | 5.467874 | 5.378189 |
| Sigma/integral-adjoint response | 5.174358 | 4.269067 | 4.221101 |
| Tile reads | 1.638344 | 0.950226 | 0.895960 |
| Input panel generation | 1.664074 | 1.564938 | 1.446401 |
| Packing itself | 0.047016 | 0.043376 | 0.041341 |

The consistent scoped gain is primarily removal of selected-lane clearing:
the first candidate already removes about 0.69 s from read callbacks. Final
read time is 45.3% below the saved parent, and the physical HVP scope 14.3%
lower. Do not attribute the full gain to the small packing scope. The extra
packing changes improve the allocation/copy dataflow and input-generation
scope, but do not independently establish a large native speedup over the
first candidate.

Counts remain **129 sweeps / 1171 batches / 1120445436 logical source
elements / 1068 product projections**. This is faster execution of the same
actions, not fewer HVPs, fewer physical contractions or a new complexity
reduction. The remaining channel compute and output scatter are still open
optimization targets.
