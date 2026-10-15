# FeCl2: inverse defects before learning a new direction

## What the diagnostic measures

Job 1249694 adds a full-audit measurement immediately after obtaining a
new direction's exact images and **before** the coarse observer absorbs it.
Let L=H-omega G, q be the admitted primal direction, and P the existing
balanced approximate inverse. Measure

$$
e=PLq-q,\qquad
\varepsilon_q=\sqrt{\frac{e^TGe}{q^TGq}}.
$$

This probes an actually encountered, previously unlearned direction. It is
not a condition-number bound. The complement fraction below uses Euclidean
projection against the existing AH basis; it is not a physical norm split.
No HVP is added. Full diagnostics add one inverse and one metric action
per expansion; scalar/default paths do not execute these operations.

The job passed all five units before chemistry and converged in 11 macros,
222 HVPs at -2181.624464858668 Eh. Its full-audit time is not a production
speed measurement, and its trajectory must not be equated to a native run.

| Macro | Unknown-direction tests | Median physical inverse error | Median Euclidean complement/error ratio |
|---|---:|---:|---:|
| 8 | 21 | 3.371 | 0.766 |
| 9 | 22 | 3.892 | 0.603 |
| 10 | 24 | 2.110 | 0.665 |
| 11 | 25 | 1.910 | 0.586 |

The new-direction inverse error is substantial even though previous
covered-space tests give near-unit spectra. This establishes a concrete
preconditioner limitation on observed new directions, not a defect in the
small Ritz solve. Small residual alignment by itself would not establish
that limitation. These errors do not exclude nonlinear keyframe curvature
variation or prove that one replacement inverse will solve the problem.

## Current-point symmetric inverse-defect correction

Let V contain the existing exact sampled primal directions, and let

$$
Y=LV,\qquad D=V-P_0Y,\qquad J=D^TY.
$$

Here P_0 is the base signed inverse at the **current** AH root, not a
previous-point secant history. For symmetric L and P_0, J is symmetric:

$$
J=V^TLV-Y^TP_0Y.
$$

When J is nonsingular, the multivariate symmetric-rank-one inverse is

$$
P_{\mathrm{SR1}}=P_0+DJ^{-1}D^T.
$$

Since D^TY=J, it reproduces the exact sampled inverse actions:

$$
P_{\mathrm{SR1}}Y=P_0Y+D=V.
$$

For numerical rank deficiency, the prototype retains the resolved signed
eigenspace of J using the existing AH dependence rule. It does not claim
full secant exactness on dropped modes. There is no new tolerance, molecule
selector, Hessian approximation in the actual AH equation, or stopping rule.

The retained balancing inverse already satisfies the same sampled inverse
identity, but its extension into the complement differs. The defect D can
contain components outside V; thus this is a different complementary-space
model, not another attempt to solve the old covered-space equations more
accurately. Neither construction guarantees improved nonlinear convergence.

## Costs and gates

The prototype uses the already retained V/HV/GV, and shared serial tall
arrays with at most the existing 60 coarse columns. It creates no full
Hessian or H/S matrix, no integral copy, and no per-thread array. Its serious
cost is applying the base inverse to every Y column at each changed root.
Consequently fewer HVPs alone cannot justify adoption: total native wall
time and memory must be checked too. No cross-root approximation/cache is
introduced to hide this cost.

Independent dense tests cover positive/negative base inverses, a nonidentity
metric and shift, exact sampled inverse actions, symmetry, and a zero-defect
exact base inverse. Job 1249700 gates those and four other unit executables
before the original-guess FeCl2 molecular trial. Controls are unchanged.
## Rejected native result

All five unit executables passed in 1249700. The original-guess chemistry
run then reached the existing 32-macro diagnostic bound without convergence:
408 HVPs, -2181.624122473344 Eh, GNU wall 18.36 s, RSS 499.33 MiB,
anode31. Application exit 2 reports this bound, not a converged endpoint.
It is inferior to archived 1249596: 12 macros / 251 HVPs / 8.30 s on the
same node, at -2181.624464871899 Eh. No old baseline was rerun.

The secant formula, symmetry and shifted dense checks passed, but those
facts do not certify a useful complementary-space inverse or improved
nonlinear trajectory. No eigenvalue cut, level shift or tolerance was
tuned to salvage this trial. The SR1 implementation, API, call site and
tests were removed. Only the unlearned-direction diagnostic remains; it
has no production/default/scalar operator overhead. No iron long run was
submitted. The goal remains incomplete.
Restored-source build-only job 1249708 passed all five unit executables
with scalar tracing enabled. It did not rerun a molecular baseline.
