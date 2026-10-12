# Next gated experiment: gradient-weighted coupled AH startup predictor

The design was implemented and tested in 1249671/1249674, then removed.
It reduced total FeCl2 actions modestly but did not improve the tail budget
exits. No production predictor API or runtime switch remains.

## Evidence and distinction from rejected trials

The retained driver exports only the last microstep as the next macro's
orbital/structure guess (`optimization/backends/ciah.cpp`). Its coarse
sample space is discarded, although its already evaluated exact samples
cover more directions than that one microstep. Full metric inverses and
small Ritz eigenproblems have been excluded as dominant tail costs.

The previous single-lowest-curvature startup seed failed (1249581).
Do not repeat it. The new question is whether the *new gradient's weights*
in the previous multi-direction model yield a useful single startup vector.
It is also not historical gradient-difference/secant fitting or copying
old HVP images into the new AH equation.

## Small predictor model

At the previous accepted-point model, retain sampled basis V and the small
symmetric curvature E=V^T H_old V. Transport its primal directions into
the new orbital/structure chart to obtain W. In the new point evaluate

$$
b=W^Tg_{\mathrm{new}},\qquad M=W^TG_{\mathrm{new}}W.
$$

Use the existing AH root-selection/dependence rules on the *approximate*
small predictor model

$$
\begin{pmatrix}0&b^T\\ b&E\end{pmatrix}
\begin{pmatrix}\alpha\\ z\end{pmatrix}
=\omega
\begin{pmatrix}1&0\\0&M\end{pmatrix}
\begin{pmatrix}\alpha\\z\end{pmatrix},
\qquad p_{\mathrm{seed}}=Wz/\alpha.
$$

E is explicitly a previous-point approximation, not the current projected
Hessian. The new workspace must evaluate H_new p_seed and G_new p_seed
fresh, exactly as it already does for its last-step startup guess. Every
subsequent Ritz equation, residual, energy/gradient keyframe and acceptance
uses the actual current-point operator. Thus a poor predictor cannot be
hidden as an inaccurate current HVP or relaxed residual certificate.

## Required gates before a molecular trial

1. Inspect the existing orbital `project_vector`/`expand_step` and coefficient
   phase-alignment APIs. Demonstrate correct primal transport, including
   coefficient normalization/phase and changed charts. Do not assume that
   a gradient pullback is a primal projection.
2. Independent dense tests must prove the small predictor formula and its
   phase/basis invariance. Test dependent transported columns with the
   existing rank rule, not a new numerical knob.
3. Verify fresh current-point images and operator counts: one startup HVP,
   no replay of every old direction. The old E must never supply a current
   AH image, norm or stopping decision.
4. Retain only one previous macro's physical orbital/structure basis and
   small curvature. Store split sectors without zero padding; never keep
   the old complete AcceptedPointContext, integrals or a per-thread copy.
   Capacity is inherited from the existing bounded coarse workspace.
5. Short original-guess FeCl2 test on Hanhai22; report reference energy as
   well as actions, macros, wall and RSS. Remove the prototype if it fails.
   No iron long run follows a failed energy or efficiency gate.

Natural-frame 1249642 is a caution: 142 versus 251 HVPs and lower runtime
looked promising, but its endpoint differed by 0.288 mEh and was rejected.
The new predictor is not allowed to redefine success around action count.

## Measured result and rejection

Build/run 1249671 used Hanhai22 anode31, 32 bound threads, exact integrals,
the unchanged analytic diagonal, chart, HVP, microstep controls and outer
criteria. All five unit executables passed before chemistry. Independent
dense tests covered generalized AH agreement, basis covariance, dependent
transported columns and fresh current-point HVP images. Restart tests
covered structure phase alignment. The first development build 1249666
failed a test because its test fixture retained a reference to a temporary
NeoProblem; it did not execute chemistry. That fixture was corrected.

| Original-guess FeCl2 | Archived 1249596 | Predictor 1249671 |
|---|---:|---:|
| Node | anode31 | anode31 |
| Accepted macros | 12 | 12 |
| Orbital HVPs | 251 | 228 |
| Final energy / Eh | -2181.624464871899 | -2181.624464855227 |
| Difference from PySCF coupled-Newton reference / Eh | +1.1834e-8 | +2.8506e-8 |
| GNU wall / s | 8.30 | 7.87 |
| GNU peak RSS / MiB | 516.23 | 491.64 |
| Complete orbital norm work / s | 0.29137 | 0.24712 |
| Complete structure norm work / s | 0.10630 | 0.10359 |
| Predictor construction / s | 0 | 0.15680 |
| Last micro-model gradient norm | 1.02985e-4 | 3.73360e-4 |
| Final physical orbital gradient norm | 1.10852e-4 | 3.69804e-4 |

The single-run time reduction is 5.2%, and the HVP reduction is 9.2%.
They are not evidence of a resolved tail: macros 9--12 still exit on
budget, and their terminal complete model residuals are larger. Both
endpoints satisfy the existing criteria and agree in energy, but this
experiment's purpose was better tail enrichment, not a less accurate
endpoint reached under the same stopping rules. No thresholds were changed.

Regression array 1249674 also converged. Archived comparisons are 1249602;
241/MnF2 nodes differ, so these times are not matched speedup claims.

| Case | Node old/new | Macros old/new | HVPs old/new | Energy difference new-old / Eh | Wall old/new / s | RSS old/new / MiB |
|---|---|---:|---:|---:|---:|---:|
| F2 | anode31/anode31 | 3/3 | 23/23 | 0 | 0.35/0.31 | 17.61/18.38 |
| 241 | anode07/anode43 | 3/3 | 21/20 | +3.8e-11 | 2.02/2.07 | 861.08/852.17 |
| MnF2 | anode07/anode43 | 10/10 | 137/138 | -3.41e-10 | 4.79/4.70 | 508.86/542.21 |

Conclusion: an approximate multivariate old-point model compressed into
one freshly checked startup direction does not cure this tail. This does
not disprove genuine current-point block enrichment or prove a specific
unrepresented eigenmode. The existing scalar evidence still places the
unresolved residual primarily in the orbital complement, not in the small
eigensolve or structure-metric inverse. No new iron task was submitted.
The experiment's source was deleted, leaving only this evidence and frozen
test artifacts; the user-owned RI edits were preserved.
Restored-source build-only job 1249679 subsequently passed all five unit
executables with scalar diagnostics enabled. It ran no molecular baseline.
