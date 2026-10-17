# FULL CI: remove the inner orbital representation

## Mathematical necessity

For a complete target-spin sector, write F=T W, G=W^T W=L L^T and
U=W L^-T. The block factors L mix only structures with the same spatial
occupation pattern. U is an isometry onto that complete spin sector.
Since the common spatial exterior map T commutes with the sector projector,
J=U^T T W is square and invertible. The selected-space CI problem becomes

    H C=M C E  <=>  (U^T H_orth U) x=x E,  x=J C.

Thus orbital exterior transforms and M actions are unnecessary **inside**
FULL CI iterations. Entry and exit maps remain necessary. A sigma uses the
existing Slater--Condon kernel between U and U^T. The existing complete-sector
inverse preconditioner transforms to U^T D^-1 U, preserving its physical
correction, root-energy shift and denominator floor. The selected-space
dimension is unchanged: 226512 here, not 853776 product amplitudes.

The common Davidson implementation now accepts an exact identity-overlap
contract. It omits the third stored M-basis, its products and its separate
reorthogonalization updates, without increasing the original subspace capacity.
The deterministic original initial block, including its random auxiliary
columns, is transformed in full. At exit, C=J^-1 x, HC=J^T h and MC=J^T x
reuse the retained Ritz images. No final sigma is needed for these images.

This is not a global dense whitening: L remains in occupation blocks and J is
never formed. Numerical carriers use the existing RHS chunk bound and shared
sigma tiles; there is no per-worker carrier copy. Admission remains the
existing bounded complete-sector factor/workspace contract. The general tiled
orthogonal-CI action is still used for nonclosed selections and larger cases
where those factors do not fit; orthogonal CI itself is not disabled.

The derivation, preconditioner congruence and exit pullback are in Section 10.1
of `article/source_column_hamiltonian.tex` and its compiled/inspected PDF.

## Qualification and stopping precision

Hanhai22 build/test job 1262779 passed 12/12 tests, including independent dense
references for nontrivial recoupling blocks and one-column RHS chunks, pure
singlet/triplet sectors, H/M exit images, metric isometry and the preconditioner
congruence. Shared Davidson tests include identity and nonidentity metrics,
one/multiple roots, recycled starts and bounded thick restarts.

No numerical tolerance, orbital coordinate, HVP or outer NEO policy changed.
The physical CI stationarity budget remains 1e-7 Eh. In isometric coordinates
it is directly the Euclidean Ritz residual, so a separate physical-residual
sigma is unnecessary. Algebraic backward-error gates use the same numerical
tolerances in the new coordinates; these norms are coordinate dependent.
Consequently the final CI solve now uses 27 actions rather than 39 and ends
at physical residual 6.34373535e-8 Eh, within the unchanged budget. The old
representation had continued to 2.06e-10 Eh. This is not identical final
state-vector precision, and the whole speedup must not be attributed solely
to faster individual actions.

## Native molecular one-step comparison

Original LOFLEA geometry and HAO orbital guess, CAS(12,12), FULL, cc-pVDZ,
libcint packed-exact integrals, NEO/Davidson, analytic orbital Hessian diagonal,
gradient 1e-3 and energy 1e-7 Eh. One accepted macro step, termination
`max_iterations` with application exit 2: **not a converged optimization**.

| Version / native job | HVP / response | Final CI actions | SCF stage (s) | Program timer (s) | Process wall (s) | Peak RSS (GiB) |
|---|---:|---:|---:|---:|---:|---:|
| saved GitHub / 1262623 | 8 / 8 | 39 | 74.965957 | 80.973591 | 93.77 | 3.373 |
| previous current / 1262735 | 8 / 8 | 39 | 100.716525 | 106.933035 | 125.88 | 3.517 |
| isometric CI / 1262780 | 8 / 8 | 27 | 78.151803 | 84.403653 | 104.20 | 3.354 |

Each uses 13 microsteps, 3 keyframes and one accepted step. All measurements
are on anode50, 8 actual OpenMP cores within one LLC, BLAS 1, 16 CPUs allocated,
48 GiB request. The new and previous-current jobs use CPUs 80--87; the saved
GitHub FULL run used 176--183. These are separate single shared-node runs,
not isolated repeat averages. New native and scalar/phase jobs ran concurrently
on disjoint CPU sets. Whole-process timing includes the roughly 158 MiB native
structure-weight output; it is not replaced by the shorter program timer.

SCF time decreases 22.4% versus the previous current implementation; process
wall decreases 17.2% and peak RSS decreases 4.6%. The new SCF time remains
4.25% above GitHub, and process wall remains 11.1% above it. **FULL has not
yet beaten the saved remote end-to-end baseline.** The remote was never rerun.

Initial energy is -398.880685140388 Eh, agreeing with the saved remote to the
printed digits. After one step, the new energy is -399.322366416769 Eh, versus
-399.322429054264 Eh remotely: +0.06264 mEh. The new accepted orbital-step norm
is 0.3599442862 versus 0.3599267076 remotely; these are different orbital points.
The endpoint difference is therefore not a same-point CI energy error, but
neither is it a converged-energy equivalence test. Physical CI accuracy and
finite optimizer trajectories are recorded separately; no claim of recovered
converged-energy parity is made from a one-step run.

The nonclosed-space regression, current-only COV job 1262787 on anode50,
retains 6 HVPs, 6 responses and 99 final CI actions. SCF/process times are
47.521203/51.91 s and RSS 2.487 GiB, versus saved GitHub 53.464025/58.06 s
and 2.612 GiB. Its after-step energy is -399.074297711215 Eh, 4.59e-8 Eh
above the saved GitHub endpoint. Thus the shared-solver changes do not show
a COV regression relative to the remote baseline in this short test.

## Current-only phase evidence

Job 1262751 qualified the preceding source with additional scalar timers:
two inverse setups total 0.307 s, two physical-residual certificates total
0.196 s, 84 inverse actions total 7.653 s. Factor setup was not the dominant
cost; merely caching it would have addressed the wrong bottleneck.

| Scope | Previous-current phase 1262736 (s) | New phase 1262781 (s) |
|---|---:|---:|
| Initial CI solve | 21.194280 | 9.477505 |
| Candidate relaxed CI solve | 16.058614 | 6.755786 |
| Trial forward, 2 calls | 39.511542 | 18.231914 |
| NEO macro | 51.130448 | 41.496798 |
| Orbital HVP, 8 calls | 8.151278 | 8.190919 |
| Coupled structure, 8 calls | 14.004684 | 13.919285 |

The two CI solves total 37.253 -> 16.233 s (2.29x), with no orbital-update or
HVP-count reduction. Parent and child timers overlap and must not be summed.
Program-wide left/right exterior-transform calls decrease from 6652/262 to
2620/101; these counters also include non-CI consumers. These are application
phase clocks, not sampled-profiler percentages. The earlier VTune ptrace-policy
failure remains the applicable sampling limitation.

Intermediate jobs 1262767/1262768 tested a version before the full initial-block
fix and identity-storage removal; they are excluded from the primary comparison.
Build 1262774 was interrupted; its zero-byte `eigen_response.cpp.o` caused the
1262775 link failure. Only that generated empty object was removed, and 1262779
rebuilt and passed all tests. This build incident is not an SCF numerical failure.

Machine-readable native measurements are in `full_ci_isometry_scf.csv`.
`full_ci_isometry_evidence.tar.gz` retains source patch, build/configuration,
tests, native/phase logs, trace, input/binary hashes, CPU affinity, accounting
and compact FULL stdout extracts. Full native stdout is kept locally and on
the cluster, and phase FULL stdout on the cluster, with hashes recorded.
