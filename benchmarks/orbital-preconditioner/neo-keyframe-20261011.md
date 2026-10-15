# PySCF keyframe study

## Scope

Study PySCF's approximate augmented-Hessian candidates and keyframe gradient
correction without importing orthogonal-orbital rotation into nonorthogonal
VBSCF. The trial prototype below was tested, then removed because it worsened
MnF2. Production retains the candidate strategy from `d86a348`; this study
does not claim a successful new solver or fewer HVPs.

Primary sources:

- [PySCF coupled Newton CASSCF](https://pyscf.org/_modules/pyscf/mcscf/newton_casscf.html)
- [PySCF one-step CASSCF](https://pyscf.org/_modules/pyscf/mcscf/mc1step.html)

These drivers have different keyframe conditions. The coupled driver's
factor-three vector-error test is not the same as the one-step driver's
conditional rollback test. Neither is simply a comparison of gradient
magnitudes. Successful keyframes correct the predicted gradient used by
subsequent microiterations; this is not equivalent to rejecting an early
candidate and demanding a fully certified solve at the original point.

## Gradient comparison

For a fixed source chart with coefficient lift $J$, our raw retraction is

$$
R(p)=x+Jp.
$$

If the raw trial gradient includes the derivatives of orbital normalization
and inactive projection, the gradient of the retracted objective is

$$
\nabla_p E(R(p))=J^T\nabla_x E(R(p)).
$$

This uses the source chart's covector projection, not tangent-coordinate
recovery and not a projection into a newly constructed trial chart. An
independent finite-difference test now checks this identity at a finite step
for partial and full orbital supports. It is a geometry chain-rule test,
not a new molecular electronic-gradient validation.

The production trial evaluator also reselects the support-preserving gauge.
The resulting raw gradient is evaluated at the gauge-transformed coefficients.
Comparing it with a source-chart prediction requires a consistent pullback;
one cannot silently treat the two raw coordinate systems as identical.

The prototype therefore omitted gauge reselection for early candidates,
compared covectors in the fixed source chart, and used the coupled PySCF
vector-error condition with our source physical dual norm. Regular certified
trials still reselected the gauge. This is a change in update strategy,
not merely a new stopping test.

## Matched native measurements

Job `281467` ran the baseline and prototype on the same node for each pair,
with identical inputs, 32 compact-affinity OpenMP cores, one BLAS thread,
exact Libcint integrals and Davidson. Outer thresholds were unchanged:
gradient 1e-3 and adjacent energy change 1e-7 Eh. Time is application
end-to-end wall time; RSS is GNU-time process peak memory. Points count
all accepted exact-gradient points, rather than grouping microiterations.

| System | Points baseline / prototype | HVP baseline / prototype | Structure actions baseline / prototype | Wall s baseline / prototype | RSS MiB baseline / prototype |
| --- | ---: | ---: | ---: | ---: | ---: |
| MnF2 | 14 / 25 | 180 / 467 | 487 / 890 | 4.317 / 9.055 | 480.8 / 493.6 |
| FeCl2 | 3 / 3 | 287 / 287 | 1007 / 1008 | 5.872 / 5.760 | 400.7 / 387.4 |
| 241 | 6 / 6 | 10 / 10 | 119 / 119 | 2.168 / 2.174 | 921.7 / 926.7 |
| 240 | 8 / 8 | 13 / 13 | 0 / 0 | 3.485 / 3.477 | 1689.5 / 1845.7 |
| F2 | 5 / 5 | 11 / 11 | 20 / 20 | 0.139 / 0.141 | 58.0 / 59.9 |

All pairs converged under the unchanged requested criteria. MnF2 energies
were -1348.893351870673 and -1348.893352949619 Eh, respectively; physical
gradient norms were 0.000611879 and 0.000842477. The energy difference is
1.079e-6 Eh: passing an adjacent-step energy criterion is not a guarantee of
agreement in final energies to 1e-7 Eh. The unchanged 240 trajectory does not
execute this coupled keyframe guard; its single-run RSS variation is not
evidence of memory allocated by the guard.

Executable SHA-256:

- baseline: `ea0909541ec29166ccb0db93a5de74985ada4040ef92e151146bb58840d2919b`
- vector prototype: `07fca4b75e09d1893a7ca19d80025ee05d0dd23b9451c5d441461ef5e540f508`

Parsed data: `results/studies/locality-20261011/keyframe-matched.csv`.

## Gauge ablation

Job `281473` retained omission of gauge reselection but used the previous
magnitude-style acceptance condition, evaluated in the source metric.
This isolates a magnitude-policy variant; it is not an exact reproduction
of the baseline's norm evaluation in the trial's newly constructed chart.

| MnF2 variant | Points | HVP | Structure actions | Wall s | RSS MiB | Final energy Eh |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Baseline, matched job 281467 | 14 | 180 | 487 | 4.317 | 480.8 | -1348.893351870673 |
| Vector guard, matched job 281467 | 25 | 467 | 890 | 9.055 | 493.6 | -1348.893352949619 |
| Magnitude guard without gauge reselection, job 281473 | 21 | 461 | 669 | 9.470 | 484.1 | -1348.893352919285 |

The magnitude/no-gauge variant also regresses. Thus the vector guard alone
does not explain the regression; the combined no-gauge/source-norm change
already changes the trajectory substantially. No parameter search was done.
The scalar-instrumented run of this ablation followed a different 25-point
trajectory (467 HVPs), so its timing is not substituted for the native row.
Its diagnostic `reliable` field records the vector test, not the temporary
magnitude test actually used for acceptance. That ablation-only discrepancy
was removed together with the prototype, not left in production diagnostics.

The vector-prototype trace rejects only the first tested early candidate:
actual norm 23.686, predicted norm 5.877, error norm 23.839. All remaining
tested early keyframes pass. A single initial rejection can change the
nonlinear trajectory; rejection count alone cannot establish its cause.

## Decision

Remove the prototype, its unused helper and helper-specific tests. Retain
the finite-step source-chart pullback regression test and include it in the
Slurm audit's five unit-test targets. Jobs `281472` and `281473` passed all
five tests on compute nodes.

After removing the prototype, job `281474` again passed all five unit tests.
Its native MnF2 run recovered 14 points, 180 HVPs and 487 structure actions,
4.583 s and 493.1 MiB peak RSS, with energy -1348.893351870673 Eh. The
restored executable hash exactly matches the immutable baseline above.
Parsed results are in `results/studies/locality-20261011/keyframe-restored.csv`.

The next mathematically justified implementation needs consistent gauge
pullback and a genuine retained-microiteration model update. It must also
distinguish a joint orbital/structure gradient from our fully relaxed trial
gradient before importing a coupled-AH prediction-error condition. Merely
adding a stricter rejection test is not the PySCF algorithm and is not
supported as a performance improvement by these measurements.
