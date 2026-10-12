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
