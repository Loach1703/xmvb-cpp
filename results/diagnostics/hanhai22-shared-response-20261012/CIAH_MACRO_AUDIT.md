# CIAH macrostep consistency diagnostic

## Purpose

Use the existing porphine triplet low-energy starting point to distinguish
the last AH Ritz residual, the accumulated frozen-model Newton residual,
keyframe corrections, and the actual fixed-C/relaxed endpoint gradients.
No HVP formula, preconditioner, coordinate map, acceptance test, or stopping
threshold is changed. No old-version run is submitted.

Hanhai22 job **1248292**, 32 CPUs, two diagnostic macrosteps. The updated
`test_ciah` passed on the compute node before the chemistry run started.
The new unit test checks accumulated/scaled H images and retained-prefix
images against a small dense reference after a gradient refresh, without
additional Hessian actions.

## Events

- `ciah_micro`: additionally records component clipping and its scale.
- `ciah_keyframe`: additionally records whether the correction was accepted.
- `ciah_micro_end`: records the actual exit reason, kept predicted gradient,
  frozen-model gradient, and accumulated keyframe correction before any
  descending-prefix restoration.
- `ciah_macro_model`: records the unshifted frozen Newton residual for the
  final accepted displacement, split into orbital/structure parts. Cached
  H images are reused; there is no added HVP.
- `ciah_macro_gradient`: an explicit nonlinear audit records the actual
  fixed-C gradient, its difference from the frozen model, and the relaxed
  orbital gradient, all projected into the same origin chart.

The nonlinear event requires existing opt-in switches
`XMVB_NEO_TRACE=1` and `XMVB_NEO_CANDIDATE_AUDIT=1`. It adds one fixed-C
keyframe per accepted macrostep. The normal relaxed gradient completion is
brought forward and cached, not computed twice. Diagnostic timing and
keyframe counts are therefore not a performance benchmark.

The coupled norm uses the existing exact orbital metric dual norm and the
structure diagonal dual norm. It is not labelled an exact inverse of the
complete structure overlap metric. The frozen-model residual is unshifted;
it is not a certificate for the shifted AH equation or nonlinear stationarity.
The existing CSV `kkt_residual_norm` still denotes the last Ritz solve and
must not be reinterpreted as the accepted macrostep residual.

Check energies and HVP counts against the already completed first two steps
of job 1248127, without rerunning the old executable. Extra scalar metric
norms and one nonlinear audit must not be presented as native runtime costs.

## Concurrent Molcas reference

Offline build **1248254** completed successfully. Calculation 1248206 then
failed on `GROUP`, which belongs to COORD/XYZ gateway input rather than the
explicit basis/atom grammar used here. The input generator now follows the
previous successful FeCl2 grammar: every atom listed explicitly, no symmetry
generators, and a mandatory NSYM=1 certificate before conversion. Corrected
Molcas calculation **1248287** was submitted; it is independent of the CIAH
diagnostic. No C++ optimizer changes were made for the Molcas correction.

## Completed two-step audit

Job 1248292 completed; native exit 2 is the deliberate two-step cap, not a
convergence flag. HVP counts remain 18/19. Endpoint energies differ from the
existing first two steps by only about 1e-11 Eh, confirming that the audit
did not materially change the recorded trajectory.

| Accepted macro | Source gradient | Frozen Newton residual | Actual fixed-C gradient | Linearization error | Relaxed orbital gradient in origin chart |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.01136091 | 0.01480096 | 0.04762050 | 0.04208605 | 0.00663483 |
| 2 | 0.00695461 | 0.00576313 | 0.01398770 | 0.01187466 | 0.00576176 |

For macro 1, the retained prefix has fixed-C energy -2244.7907679919244 Eh,
above the origin -2244.7907841974 Eh. The final structure reoptimization
rescues it to -2244.7907984105095 Eh. The same prefix was a previously
accepted keyframe: mismatch / predicted-gradient norm is approximately 2.84,
below the current inherited trust factor 3. Thus the code accepts a keyframe
with a much larger true joint gradient, while continuing with origin Hessian
images. The final accumulated direction loses frozen-model descent, requiring
the prefix repair. Macro 1 exits on micro budget; macro 2 exits on keyframe
mismatch and also restores a prefix.

These observations locate a failure of the finite coupled-step/model
management: reliable energy decrease after CI relaxation is not evidence
that the candidate was a good coupled Newton correction. They do not isolate
every underlying cause or prove all HVP directions correct; earlier
accepted-point derivative checks cover sampled directions only. No threshold
or optimizer policy was changed as a result of this diagnosis.

## Full comparison, not restarted step counts

Array **1248419**, submitted on Hanhai22 with 32 CPUs per task:

- `_0`: current NEO, analytic-diagonal preconditioner, from the unchanged
  `inputs/triplet/input.xmi` ($GUS initial orbitals), not an optimized endpoint.
- `_1`: PySCF default AH from the same raw $GUS orbitals, converted to a
  core/active orthonormal frame without changing those subspaces. Each program
  solves its own initial triplet CI problem. No endpoint CI vector is reused.

CAS(8,11), Cartesian cc-pVDZ, no ECP/symmetry/RI, 1e-7 Eh and 1e-3 gradients,
1000 macro safety limit, and the partition's 15-day scheduler limit are kept.
The C++ executable is the already tested diagnostic build, with scalar events
only and extra nonlinear auditing disabled. Original input/basis hashes are
saved for both tasks. Both count macrosteps from original iteration zero.
No new SCF guess is generated. This is an unchanged-policy full comparison,
not a claimed convergence improvement.
