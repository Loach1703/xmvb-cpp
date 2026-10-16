# C14 metric-CG diagnosis and bounded factor reuse

## Failure attribution

Original run `1264455` failed before accepting an orbital step. A bounded
one-macro diagnostic `1264468` identified the exact caller as the CI transport
constructor, not an orbital Hessian preconditioner or CI eigensolver failure.
The constructor solves the selected-structure Gram equation for its residual
weighted transport adjoint. In that capture the right-hand-side norm was
`2.7324861941668125e-9`; the unchanged true-residual target was
`4.0717217239242978e-17`. Jacobi PCG stopped after 858 iterations with a true
residual `3.2656144719694913e-11`, about 1.2% relative. This was not a failure
caused by demanding accuracy below the attainable floating-point level.

Independent accepted-point diagnostic `1264469` measured:

| Quantity | Value |
|---|---:|
| Selected structure dimension | 429 |
| Orthogonal direct-CI carrier | yes |
| Resident bounded selected Gram | yes |
| Relative Gram asymmetry | 7.02e-16 |
| Minimum / maximum eigenvalue | 5.1978e-4 / 7548.7173 |
| Gram condition number | 1.4523e7 |
| Diagonally scaled condition number | 1.5140e7 |
| Relative approximate-diagonal error | 0.3699 |

The metric is positive definite but ill-conditioned in the selected structure
coordinates. The diagonal supplied by `StructureAction` is an approximate
positive preconditioning diagonal, not the exact selected Gram diagonal;
its scaling does not improve this conditioning. The CG implementation also
cannot infer singularity from a cycle's increasing Euclidean residual: that
residual is not guaranteed to decrease monotonically. The observed guard is
a solver failure report, not evidence that the physical metric has a nullspace.
The cure used here is an already available exact factor, not a modified
stagnation threshold or an enlarged iteration budget.

For the independently reconstructed actual transport RHS, original PCG failed
after 1716 iterations; an LLT solve certified against the **original** Gram
action gave:

- RHS norm `2.7355527853370451e-9`, target `4.0762913008586244e-17`.
- LLT true residual `2.7926510442299701e-22`, relative `1.0208726584253757e-13`.
- An independent deterministic RHS also certified, relative residual
  `3.6072617572171043e-12`; its original PCG solve failed.

The last digits and iteration counts differ between separately reconstructed
accepted points. The failure attribution does not depend on treating their
tiny CI residual vectors as bitwise equal. Both independently expose the same
poorly preconditioned Gram solve and pass the original tolerance with LLT.

## Fix and invariants

`SelectedMetric` owns the already resident Gram and its positive LLT factor.
`PreparedMetric`, Davidson, the accepted-point context, and CI transport share
one immutable instance, rather than destroying the Davidson factor and
re-solving with only the approximate diagonal. A trial orbital chart obtains
its own factor and invalidates the cached trial whenever the chart changes.
The structure tangent uses the same exact constrained horizontal inverse
formula as the existing complete-space inverse, tested against an independent
dense restricted metric for one and two roots. Exact factors bypass redundant
Krylov recycling; the true-residual certification itself is unchanged.

The existing selected-Gram workspace contract is unchanged: only a Gram that
already fits its 2 MiB quarter-workspace is resident. For C14 the matrix is
1,472,328 bytes, and its LLT has the same order of storage, shared rather than
replicated per consumer or worker. Large selected spaces keep streamed actions
and do not acquire an unbounded selected-space or spin-string-product matrix.
No arbitrary full-CAS closure is assumed.

Orbital coordinates, HVP formulas, CIAH policy, all input settings, and inner
and outer accuracy thresholds are unchanged. Error messages now retain caller,
dimension, residual target, and iteration counts so the next failure is not
reported as an unexplained generic stagnation.

## Verification

New unit checks cover shared ownership/lifetime, an ill-conditioned selected
Gram at the original true-residual tolerance, changed trial charts, one- and
two-state horizontal inverses, and finite CI transport derivatives. Build,
tests and the full native C14 rerun use Hanhai22 job **1264476**
with 32 physical cores. Compilation succeeded; **7/7 tests passed**, 0.32 s.
All eleven changed production/test source files matched the compute-node
source manifest. The final run uses the exact same `C14.xmi` and preserves
its `ITMAX=200`; there is no added diagnostic macro-step cap.

The native rerun has reached the initial energy `-534.3944115253 Eh` and is
running its first CIAH macro step on `anode40`. Full C14 convergence has **not
yet been established**. Numerical diagnostics are evidence for this solve fix,
not a completed optimizer run. Build/test and diagnostic evidence is retained
in `metric-fix-evidence.tar.gz`; native long-run output remains in
`/home/guqqgroup/xiatao/xmvb-runs/c14-ciah-20261016/fix-1264476`.
