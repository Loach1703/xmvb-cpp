# Exact support-amplitude maps and single-pass H-image consumers

2026-10-17. Parent `b1f38c56`. No historical executable or baseline job was
rerun. Production thresholds, CIAH policy, spin space, structure space and
exact integral representation are unchanged. The unrelated local RI edits
were not transferred or edited.

## Implementation

`SpinFactorSource` now describes both its source factors and an optional
exact narrow-panel amplitude multiplication. The VB adapter implements
`Y V = T_alpha [Z (T_beta^T V)]` and its transpose with existing original
support incidence lists. No dense Z, transformed product amplitude, all-source
Gram or per-thread carrier is introduced. The raw parity contributions are
combined before the exterior transform; both signs and halves are retained.

`SpinAmplitudeWorkspace` selects between exact factor moments and support-map
multiplication using their operation counts. Integral and oriented-generator
adjoints share this implementation. The source/channel planner includes
callback scratch and accounts for the changed moment work. This is a choice
of contraction order, not an approximation, orbital-closure assumption,
physical-space admission test, or compatibility fallback.

`project_source_adjoint` executes exactly one H-image sweep. A cheap resident
right-covector accumulation is retained when it fits the same bounded budget
and avoids more generator-link work; otherwise image tiles are consumed via
the narrow amplitude products. Right-source streaming no longer encloses the
Hamiltonian sweep. Both schedules preserve the orientation of D((H-E)X,Y).
No accepted-point image cache was added. The shared numerical budget remains
64 MiB; immutable input metadata and required output arrays are separate.

Theory and complexity are recorded in Section 11 of
`article/source_column_hamiltonian.tex` and the rebuilt, inspected PDF.
The source-pair term can lose its explicit left-source factor for sparse
support, but the independent dense channel-Gram term is not eliminated.

## Qualification and independently checked source path

Final qualification **1265522**, anode62: **14/14 CTests passed**, 6.37 s.
The broadened three-structure test compares source and product orders for
arbitrary nonclosed support, unequal spin spaces and nonstationary responses.
The tests also cover positive/negative parity, integral derivatives, HVPs,
oriented generator pairing and constrained workspaces. Callback-only tests
exercise the map without an extra source-family allocation.

Final 32-core diagnostic **1265525**, anode04, CAS(12,12), one structure:

| Independently checked quantity | Absolute error |
|---|---:|
| Active-integral gradient | 6.07e-15 |
| Directional H/M images | 1.95e-18 |
| Integral-adjoint HVP | 9.54e-18 |
| Independent Lagrangian finite difference | 6.49e-11 |

The source path is actually used. Its logs show zero regenerated left-source
columns in mapped integral/generator contractions and one H-image sweep per
image consumer. Accepted/directional product-array sizes are zero, the
direction workspace bound is 36,588,672 bytes and peak process RSS is
74,236 KiB. This is a short operator qualification, not an SCF wall-time
speedup claim. The diagnostic phase instrumentation is not used in the native
timings below. Tiny complete carriers occur only in independent test oracles.

## Requested end-to-end tests: no demonstrated broad-space speedup

Frozen bundle **1265523**, native array **1265524_[1,2,6]**. All three run to
convergence from their original input guesses on Hanhai22, 32 physical cores,
OpenMP cores/close, BLAS1, 64 GiB requested, exact Libcint/cc-pVDZ, HAO,
ISCF=9/CIAH and structure-space Davidson. No added iteration/time cap is used.
Input hashes are identical to the saved controls. COV also matches the archived
workspace-lease input byte for byte. All nodes use AMD EPYC 9654 processors.

Application wall includes initialization and SCF; GNU process wall additionally
includes final output/teardown. RSS is GNU time's process peak, not summed Slurm
RSS. Native runs contain no profiler collection.

| Space | Structures | Steps old/new | HVP old/new | Application s old/new | Process s old/new | RSS GiB old/new |
|---|---:|---:|---:|---:|---:|---:|
| COV | 132 | 4 / 4 | 36 / 36 | 17.208 / 19.475 | 17.40 / 19.72 | 2.572 / 2.752 |
| ION(0:2) | 47,256 | 4 / 4 | 44 / 44 | 104.286 / 120.715 | 107.11 / 124.01 | 2.824 / 3.031 |
| FULL | 226,512 | 4 / 4 | 40 / 40 | 102.727 / 100.683 | 112.47 / 111.18 | 3.573 / 3.578 |

| Space | New total energy / Eh | New minus saved / Eh |
|---|---:|---:|
| COV | -399.106646267329 | +4.04e-12 |
| ION(0:2) | -399.397232509244 | +9.66e-13 |
| FULL | -399.417495587054 | +9.66e-13 |

COV control is the excellent saved **1264892**, not an older slower 24-s
version. ION control is retained **1265266_1**, not one of the discarded
inverse-model experiments. FULL control is saved **1265114_2**. COV compares
anode40/0-31 to anode23/160-191, ION anode65/96-127 to anode04/160-191,
and FULL anode42/64-95 to anode23/128-159. These are single-run cross-node
comparisons: the slower COV/ION observations are retained, not hidden, and
neither their regression nor FULL's small improvement is attributed to the
new source kernels.

These broad-support examples continue to choose the established product-tile
orders. No change was made to force an otherwise expensive source order.
Thus the implemented algebra/low-memory image handling is qualified, but
**the requested broad-space end-to-end acceleration is not demonstrated**.
Further work must address the executed product-tile chain, not claim that a
source-kernel improvement automatically accelerates these examples.

## Bounded attribution, not another convergence benchmark

**1265545**, anode04, profiles only the first unchanged current ION macrostep
with gprofng plus phase instrumentation. Reaching that requested one-step
limit returns 2; it is intentional, not a failed convergence run. Native
converged array runs above all return 0. The capture reports no collection
errors/warnings. Its percentages include synchronization and unidentified
libc frames and are not presented as precise useful-compute fractions.

There are no factor-adjoint/source-image/source-sweep trace records in this
capture. The phase `ci_integral_adjoint` is zero, while the executed chain is
`contract_product_tiles` / `product_jet` / `contract_tiled`. First-step phase
walls include CI sigma 8.961 s, product generation 4.683 s, product projection
4.104 s, and left/right exterior transforms 1.953/1.708 s. These are nested,
instrumented values, not additive native time components. They identify the
remaining path; they do not explain a cross-node timing difference causally.

## Evidence

`comparison.csv` is parsed from saved CSV controls and current outputs;
`native.csv` contains only the three requested completed runs. The compressed
archive retains final source/build manifests, qualification logs, frozen
inputs, complete native outputs and traces, the source diagnostic and profiler
reports, excluding binaries and the large capture directory. The original
gprofng capture remains at
`/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015/c12-ionic-profile-1265545/capture.er`.
Initial implementation trials are separate Slurm jobs 1265493/1265507,
1265510, 1265514, 1265516/1265519 and 1265521. They are not used as final
performance results. Pending fixed-node job 1265509 was cancelled before
execution; no historical job was run or cancelled.
