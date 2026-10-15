# Equal-weight multistate CIAH

## What caused the old route

`neo.cpp` selected CIAH only for one root. `ciah.cpp` and `CiTransport`
explicitly rejected multiple roots. State-averaged calculations therefore
used `ResponseNeoWorkspace`, not the current keyframe algorithm. This was
an implementation coverage gap, not a mathematical need for a second solver.
The old branch unconditionally required lowest-curvature certification and
evaluated the full physical orbital dual norm during candidate refinement.
It prioritized `structure_curvature` expansions when that certificate failed.

Saved array 1255187: OEO step 1 used 42 orbital HVP and 514 structure actions;
step 2 used 23 and 871, including four rejected trials. BLO step 1 used 39 and
1097, also with four rejected trials. These are action counts, not profiler
time percentages. Do not loosen tolerances to suppress them.

## State-subspace objective

For r equally weighted states, let C be the raw coefficient block, M the
structure overlap, H the structure Hamiltonian, G=C^T M C and K=C^T H C.
The fixed-subspace electronic objective is

\[
E(C)=r^{-1}\operatorname{Tr}(G^{-1}K),\qquad
\nabla_C E=2r^{-1}(HC-MC G^{-1}K)G^{-1}.
\]

Independent column normalization does not implement this objective once
finite-step states cease to be mutually M-orthogonal. The new keyframe uses
an r-by-r Cholesky factor and projected Hamiltonian eigensolve. With T such
that T^T G T=I and T^T K T=diag(e), use CT in the existing orbital adjoint,
then pull the coefficient derivative back as

\[
\nabla_C E=2r^{-1}[HCT-MCT\operatorname{diag}(e)]T^T.
\]

No full structure eigensolve or dense structure metric is added. This equals
the Rayleigh quotient at r=1 and is invariant under any nonsingular change
of the equal-weight state frame. Unequal weights are outside this extension.

## Transport and pullback

The existing selected-space finite projection is applied to the entire
state block. The cross-Gram right-hand side and metric solves are shared
where their block interfaces support it. Adjoint contractions are summed
over states. This is a projection into the supplied VB structure space;
closure of that space under orbital transformations is never assumed.
The structure coordinate scale is sqrt(2/r), consistently in the origin
gradient, finite-gradient pullback and residual-weighted Hessian correction.

On restarting from relaxed CI roots, a small orthogonal Procrustes rotation
aligns the recycled structure direction with the new state frame. At r=1
this reduces to the former sign correction.

## Validation contract

- Finite-difference orbital and coefficient gradients for full and arbitrary
  selected structures, including nonorthogonal supplied CI columns.
- Invariance under nonsingular state mixing and uniform rescaling.
- Finite selected-space transport gradient and residual-weighted HVP
  correction with two states, independently differentiated.
- Existing single-state transport, CIAH, metric and Hessian tests.
- Small state-average molecular tests before rerunning pyrrole OEO/BLO.

The 7.0103 mEh historical OEO difference remains a separate unresolved
question: block-LBFGS does not use these CIAH keyframes. Historical independent
PySCF fixed-orbital results reproduce the old lower branch with spin squared
0.75 for both states. The new attachment has a different initial energy
(-417.496351311941 versus the historical -417.492571494774 Eh).
Job 1255282 failed at import because the diagnostic script shadowed the
PySCF package; the renamed script in 1255293 completed successfully.
From the new attachment its initial average is -417.49635131423315 Eh
(2.29e-9 Eh from CPP); final average is -417.49785390489683 Eh, with roots
-417.5152841118022 and -417.4804236979914 Eh. Both spin squares are 0.75.
The optimizer took four macro iterations / 58.68 s. Thus differing inputs
do not explain the remaining low CPP branch; its cause is not yet identified.

Build 1255278 passed 13 unit executables including multistate transport and
keyframe finite differences. Final build 1255305 also passed all 13 tests,
after removing the retired response solver and testing state-frame alignment.
1255297 was a mirror cleanup/configuration failure; 1255299 caught a test-only
matrix/vector type mismatch, corrected before the successful final build.
New arrays 1255306 (two small SA regressions and coupled derivative checks)
and 1255310 (pyrrole OEO/BLO NEO) use frozen build 1255305. At submission no
end-to-end speedup or molecular convergence has been verified.

Array 1255306 completed both small regressions. F2-SA-HAO converged in
4 macros (0.413738 s native wall), F2-SA-OEO in 2 (0.180319 s).
The independent coupled orbital/structure HVP finite-difference relative
errors were 4.96e-9 / 1.51e-9 (HAO), 7.14e-9 / 4.12e-9 (OEO).
The pyrrole jobs 1255310_2/3 are running; their total benefit remains pending.

## First molecular results (1255310)

BLO converged in 2 macros / 13 orbital HVP / 13 structure actions,
141.866773 s internal wall, 144.77 s GNU wall, 24099484 KiB peak RSS.
Its state energies are -417.513869168244 and -417.478838172582 Eh;
average -417.496353670413 Eh. This is slower than block-LBFGS (63.37 s),
despite far fewer actions than the old NEO route.

First-macro orbital HVP / structure actions:

| Case | Old coupled solver | Unified CIAH |
|---|---:|---:|
| OEO | 42 / 514 | 18 / 18 |
| BLO | 39 / 1097 | 10 / 10 |

New first macros have zero rejected trials. OEO macro 1 reaches
-417.497853900970 Eh (printed -417.4978539010), near the same-guess PySCF
terminal energy. Macro 2 decreases further to -417.4978572749 Eh, with
step-Hessian-step -3.3430639744e-6 Eh and actual reduction 3.3739095215e-6 Eh.
Its physical gradient grows from 8.61e-5 to 6.98e-3. This is a downhill
negative-curvature step, not evidence that the new run has converged to the
PySCF endpoint. No claim of a resolved 7 mEh difference is justified yet.
At this snapshot OEO has 2 accepted macros / 24 HVP and is still running.
