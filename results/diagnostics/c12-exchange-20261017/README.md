# C12 physical response: exchange reduction and support reuse

## Mathematical changes

1. For complete left/right carriers of the same spin-exchange parity in
   identical spin-string spaces, the beta orbital-generator cotangent equals
   its alpha counterpart. Sum twice the alpha contraction instead of gathering
   and contracting beta spectators.
2. Tiled sigma emits weighted representatives, not automatically parity
   carriers. Complete resident output panels are first projected with
   `(I+sP)/2` in place. Each unordered spin product has one owner. Partial,
   evicted output panels retain the general two-spin contraction. No carrier
   buffer, larger workspace, precision cutoff, or physical-response omission
   is introduced.
3. Directional tiles borrow the accepted immutable compact CSR support when
   panel extents agree. They do not regenerate identical symbolic support.

Derivation and applicability: [exchange-reduced generator](../../../article/exchange_reduced_generator.md).
These changes preserve the physical orbital HVP, CI Hessian/metric and both
directions of orbital-CI coupling, as well as the CIAH policy and tolerances.
They apply equally to exact/RI active integral data; this run tests exact.

## Qualification

Hanhai22 Slurm **1264598**, `anode40`, 32 physical cores, affinity `60-91`,
GCC 14.1.0, Release C++23, `-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`,
BLAS threads=1. Job completed exit 0. All five modified production/test source
files match the recorded source hashes. Unrelated local RI edits were not
transferred. Existing phase timers attribute scopes; this is not a sampling
profiler capture, and nested scope times must not be summed as disjoint shares.

**9/9 CTests passed**, 3.42 s. New independent generator tests cover both
exchange parities, incomplete left panels, multiple right arguments,
transposes and symmetric pullbacks. A throwing beta spectator reader confirms
those reads are absent. Existing structure tests cover nonclosed selected
spaces, two states, unequal spins and narrow evicting workspaces.

An initial candidate (job 1264595) incorrectly treated weighted sigma images
as parity carriers. The independent structure-response test rejected it;
no SCF timing from that candidate is used. The in-place image projection
above fixes that exact issue without changing any test tolerance.

CAS(12,12), 132-structure independent checks:

| Check | Error norm |
|---|---:|
| Selected H action | 1.97942125e-14 |
| Selected M action | 1.90486032e-15 |
| Structure-response image | 3.36632466e-16 |
| Structure-response integral adjoint | 1.49791151e-15 |
| Product-jet HVP vs gradient finite difference | 5.91078001e-11 |
| Lagrangian directional derivative | 5.21188642e-12 |

`product-jet12` reports `reference_checked=0` for its optional complete-carrier
replay; its independent finite-difference checks still execute and pass. Zero
placeholder replay errors are not numerical evidence.

## Native end-to-end result

Same byte-identical original LOFLEA input/guess, CAS(12,12), COV 132, HAO,
exact cc-pVDZ, ISCF=9 and Davidson. "C12" denotes the CAS12 fixture and
this matching molecule, not a separate user-provided C12 geometry.
Baseline **1264591**, a652a728, is saved data: no old executable was rerun.

| Metric | Saved baseline | New |
|---|---:|---:|
| Initial energy / Eh | -398.533053417516 | -398.533053417517 |
| Final energy / Eh | -399.106646267333 | -399.106646267330 |
| Accepted steps | 4 | 4 |
| Physical orbital HVPs | 36 | 36 |
| Physical structure responses | 36 | 36 |
| SCF wall / s | 33.136628 | 32.038107 |
| Application end-to-end / s | 36.085948 | 34.899360 |
| Process peak RSS / MiB | 2692.195 | 2697.812 |
| Mean accounted CPU cores | 27.67 | 27.67 |

Native wall time decreases **3.29%** (1.034x); RSS increases 5.62 MiB (0.21%),
so memory is effectively unchanged, not improved by a new scaling law.
Both terminate with `ciah_projected_gradient_tolerance`. Finite energy
differences are about 3e-12 Eh; the changes do not reduce the HVP count.

## Separate phase-enabled runs

| Inclusive scope | Saved | New |
|---|---:|---:|
| 36 coupled actions / s | 16.467049 | 15.155833 |
| Physical orbital HVP / s | 9.853362 | 9.266741 |
| Structure-side scope / s | 6.593151 | 5.867903 |
| Generator adjoint: 704 calls / s | 1.075856 | 0.384415 |
| Symbolic support construction | 812 / 1.324120 s | 524 / 0.419583 s |
| Product projection / s | 1.682614 | 1.800119 |
| Product generation / s | 2.487389 | 1.500037 |
| 88 physical sigma/adjoint sweeps / s | 14.789951 | 13.908057 |

The in-place parity projection adds work to the projection scope; generator
and support savings exceed that cost. Both phase runs still have 36 HVPs and
the same 165 tiled sweeps / 3109 batches / 1,964,120,314 logical source elements.
Thus this is an exact contraction reduction, not a claim that channel GEMM or
sigma halo asymptotic scaling has changed. Remaining tile compute time is
9.45 s in this diagnostic run (9.28 s saved), making it the next larger target.

Raw logs, input/binary/source hashes, CTest output and build flags are retained
in `run-1264598/`; `comparison.csv` and `phase-summary.json` are parsed data.
