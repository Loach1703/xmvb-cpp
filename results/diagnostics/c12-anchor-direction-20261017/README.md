# Accepted carrier: exact orbital-direction reuse

Theory: [accepted carrier direction](../../../article/accepted_carrier_direction.md).
The retained exterior carrier generates its tangent through the one-body
generator K=dot(R) R^-1. A triangular solve supplies K. There is no explicit
inverse, complete spin-square generator, larger cache, per-worker coefficient
buffer, or dropped physical response. Same-parity inputs use alpha links and
an in-place exchange projection; unequal spins retain both contributions.

Hanhai22 job **1264605**, `anode40`, 32 physical cores, same byte-identical
LOFLEA CAS(12,12), COV132, HAO, exact cc-pVDZ, ISCF=9 and Davidson. Existing
saved baseline 1264598 is not rerun. Unrelated local RI changes are not synced.
GCC14.1, Release C++23 and BLAS threads=1 are unchanged.

**9/9 regression tests passed** (3.48s). New panel tests compare general and
alpha-only generator actions against independent full-carrier actions on equal
and unequal spin spaces, sparse and dense immutable columns and strided outputs.
Physical C12 HVP finite-difference error: **5.91078001361e-11**; Lagrangian
directional-derivative error: **5.21188642374e-12**.

| Native metric | Saved 1264598 | New 1264605 |
|---|---:|---:|
| Accepted steps | 4 | 4 |
| Orbital HVP / physical CI responses | 36 / 36 | 36 / 36 |
| Final energy / Eh | -399.106646267330 | -399.106646267331 |
| SCF wall / s | 32.038107 | 31.425332 |
| Application wall / s | 34.899360 | 34.213118 |
| Process peak RSS / KiB | 2762560 | 2727396 |

The single native run is 1.97% faster, with 34.34 MiB lower process peak RSS.
This is not evidence for a large kernel speedup: the separate phase-enabled
runs are 35.218865s saved versus 36.111756s new. Their 36 coupled actions cost
15.155833s versus 15.721623s, while the physical sigma/adjoint scopes are nearly
unchanged (13.908057s versus 13.920101s). Directional exterior calls decrease
from 904 to 760 and are replaced by 144 bounded generator calls (0.172436s).
The new input projection and generator have costs of their own. No claim of
consistent phase-wide acceleration is made from these measurements.

The unchanged 165 tiled sweeps / 3109 batches confirm that the Hamiltonian
channel multiplication is not eliminated by this reuse. Next investigate its
integral Gram multiplication: current row-stripe workers repeatedly read the
whole right-channel panel. Two-dimensional output tiling can reduce those
repeated reads without changing the Gram, adding numerical buffers, or tuning
a physical tolerance. Existing versions must not be rerun.
