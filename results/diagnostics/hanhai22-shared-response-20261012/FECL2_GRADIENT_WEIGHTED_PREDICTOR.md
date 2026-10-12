# Next gated experiment: gradient-weighted coupled AH startup predictor

This is a design, not an implemented or benchmarked improvement.

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
