# Necessary maps before faster maps

## Mathematical boundary

The finite CI transport projects a target wave function into the selected
VB space at the displaced orbitals. For a nonclosed selection, its restricted
normal equation and adjoint pullback cannot simply be deleted. They are
necessary for this particular finite coordinate map, not a universal
requirement of VBSCF. Two consumers at exactly the same macro coordinates,
however, need only one evaluation of that immutable map.

In a complete spin sector the common spatial exterior transformation commutes
with the sector projector. Its Moore--Penrose identity cancels the interior
forward/inverse maps in the carrier-diagonal preconditioner. The selected-state
horizontal constraints then require only two state-rank projections using the
already available M C, rather than two coordinate-metric CG solves. Neither
the projection, Riesz maps, root-energy shifts nor signed denominator floor
is discarded. No complete-space cancellation is used for COV or an arbitrary
nonclosed selected space.

The derivation is Section 10 of `article/source_column_hamiltonian.tex` and
its compiled PDF. The obsolete composite carrier-preconditioner entry point
was removed, not retained as an alternate implementation. Existing workspace
and factor bounds are unchanged; there are no per-thread carrier copies.

## Qualification

Hanhai22 compute job 1262723 built the changed source and passed 11/11 tests.
Independent dense references check the complete-sector inverse, both horizontal
Riesz maps, single/equal-weight two-state constraints, and symmetry. Additional
tests check complete singlet and triplet sectors that do not span the whole
fixed-Ms product carrier, and multiple RHS chunks. Transport and keyframe
derivative regressions also pass. No threshold or optimizer policy changed.

## Native one-step molecular measurements

Same original LOFLEA CAS(12,12) HAO input, cc-pVDZ, libcint packed-exact
integrals, ISCF=9, Davidson, analytic Hessian diagonal, 1e-3 gradient and
1e-7 Eh energy tolerances. Exactly one accepted macro step is requested;
all four rows end with `max_iterations` / exit 2, not convergence.

| Space | Version / job | HVP / response | SCF stage (s) | Process wall (s) | Peak RSS (GiB) |
|---|---|---:|---:|---:|---:|
| COV, 132 structures | saved GitHub / 1262622 | 6 / 6 | 53.464025 | 58.06 | 2.612 |
| COV | cancelled maps / 1262709 | 6 / 6 | 46.721394 | 51.21 | 2.512 |
| FULL, 226512 structures | saved GitHub / 1262623 | 8 / 8 | 74.965957 | 93.77 | 3.373 |
| FULL | cancelled maps / 1262735 | 8 / 8 | 100.716525 | 125.88 | 3.517 |

Only current executables were run in this round. GitHub 946934570c096ba5072923effa7fb063736731f1
uses the already saved measurement, never another baseline run. All jobs are
on anode50, 8 actual OpenMP cores in one LLC, BLAS 1, 16 allocated CPUs,
48 GiB request. The new native FULL job uses CPUs 80--87; the saved baseline
FULL job used 176--183. Each is a single shared-node measurement, not an
isolated repeat average. Actual affinity and native GNU-time records are kept.
COV job 1262709 preceded the last FULL-only horizontal cancellation; its
COV code path is the same as the final tested source.

Initial energies remain -398.533053417516 Eh (COV) and -398.880685140388 Eh
(FULL). New after-step energies are -399.074297734368 and -399.322429044560 Eh.
They differ from saved GitHub endpoints by 2.27e-8 and 9.70e-9 Eh, respectively,
below the unchanged 1e-7 Eh criterion. These are finite-step comparisons,
not converged-energy claims or evidence of bitwise-identical trajectories.
Final CI actions remain 99 (COV) and 39 (FULL); microsteps/keyframes remain
10/3 and 13/3. Compact machine-readable results are `cancelled_maps_scf.csv`.

COV SCF-stage time is 12.6% below GitHub. FULL is 10.7% below the previous
current measurement (112.829673 s) but still 34.3% above GitHub. The goal of
beating GitHub end-to-end is therefore **not yet achieved for FULL**. Full
stdout prints roughly 158 MiB of structure weights; process wall includes
this output and is retained separately from the application timer.

## Attribution, not native timing

Current-only phase runs 1262646 (before cancellation) and 1262736 (after)
record the following nested scopes. Do not sum parent and child timers.

| Scope | Before (s) | After (s) | Calls, both |
|---|---:|---:|---:|
| NEO preconditioner | 5.547705 | 0.747204 | 7 |
| CIAH macro | 61.856896 | 51.130448 | 1 |
| Trial forward | 49.227146 | 39.511542 | 2 |
| Orbital HVP | 8.403044 | 8.151278 | 8 |
| Coupled structure | 14.554378 | 14.004684 | 8 |

Existing scalar events in 1262736 further report initial CI solve 21.194280 s
and candidate relaxed CI solve 16.058614 s; two fixed-coefficient keyframes
have no CI solve. Orbital preparation/integral-transform fields are each
below 1 s. The CI-solve timer also includes complete-sector inverse setup;
it is not exclusively sigma time. Separating unavoidable projection from
repeated coordinate-map work in this scope is the next mathematical audit.
The sevenfold local preconditioner improvement does not establish a similar
end-to-end improvement.

VTune job 1262641 failed before executing the application because compute-node
ptrace policy blocks software collection. No sampled CPU attribution is
claimed; the timings above are application phase clocks, not profiler samples.

`cancelled_maps_evidence.tar.gz` preserves the tests, configuration, input/binary
hashes, native/phase records, affinity, source patch, and compact FULL output
extracts. The new FULL raw stdout is retained locally and on the cluster;
the phase FULL raw stdout remains on the cluster with its hash recorded.
