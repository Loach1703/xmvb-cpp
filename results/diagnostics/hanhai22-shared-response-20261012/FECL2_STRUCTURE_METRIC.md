# FeCl2: structure metric and keyframe decisions

## Definition and PySCF correspondence

Let M be the VB structure overlap, Q the origin horizontal tangent basis,
and g_s the structure covector in that basis. The structure block of the
generalized AH metric is

$$
M_s=Q^T M Q.
$$

Together with the existing orbital metric G_o, the consistent dual norm is

$$
\|g\|_*^2=g_o^T G_o^{-1}g_o+g_s^T M_s^{-1}g_s.
$$

The previous structure norm used the ambient Jacobi approximation

$$
g_s^T Q^T[\operatorname{diag}(M)]^{-1}Q g_s
$$

instead of the second term. This is a useful positive preconditioner pairing,
but it is not the dual norm induced by the complete AH structure metric.

PySCF `newton_casscf.update_orb_ci` uses Euclidean norms for its CI coordinates
in the orthonormal carrier. A nonorthogonal structure-coordinate transformation
must transform the metric as well as the gradient. For a coordinate map J,

$$
g_s'=J^Tg_s,\qquad M_s'=J^T M_sJ,
\qquad (g_s')^T(M_s')^{-1}g_s'=g_s^TM_s^{-1}g_s.
$$

Thus copying a Euclidean CI norm into VB coordinates, or replacing the complete
inverse by the ambient Jacobi pairing, does not reproduce the chosen AH metric.
This correspondence does not copy PySCF orbital rotations or assert identical
finite-coordinate trajectories. The metric here remains the immutable origin
model, consistent with the origin-chart gradient pullback and cached H images.

## Counterfactual diagnostic without changing the trajectory

Hanhai22 job **1248690**, FeCl2 original guess, 32 CPUs, exact integrals,
analytic diagonal preconditioner, four macros. The production decisions still
used the diagonal approximation. At every actual keyframe, a bounded diagnostic
assembled the 419-dimensional horizontal metric (1,404,488 bytes per matrix),
certified symmetry and positive Cholesky factorization, and recorded both norms.
No global structure H was assembled and no matrix was replicated per thread.
The dense diagnostic is gated by the existing explicit candidate-audit option
and a 16 MiB matrix admission bound; it is not a production metric path.

Eight keyframes were evaluated. Seven decisions agree. The final keyframe of
macro 4 is rejected with the diagonal norm but accepted with the full norm:

| Quantity | Diagonal approximation | Full metric |
|---|---:|---:|
| Predicted joint gradient norm | 0.0088117704 | 0.0118104367 |
| Gradient-mismatch norm | 0.0295041609 | 0.0285108545 |
| Mismatch / predicted norm | 3.3483 | 2.4140 |
| Same trust factor 3: decision | reject | accept |

The near-stationarity alternative in the keyframe test is not responsible for
this flip. None of the eight sampled gradient-decay trigger comparisons flip;
this does not certify all intervening microstep triggers or later trajectories.

## Matrix-free repair

`StructureTangentOperator::dual_norm` solves the horizontal Riesz equation using
matrix-free PCG and a positive Jacobi preconditioner. The recursive residual
does not certify acceptance: the true residual is recomputed before return.
The relative residual accuracy is sqrt(machine epsilon), the same numerical
contract already used by the orbital metric solve. A complete CG cycle restarts;
stagnation above that accuracy throws instead of returning an approximate norm.
Workspace is O(n_str) vectors, not O(n_str^2), and no thread owns a metric copy.

The CIAH workspace, microsteps, and keyframe tests now use this full combined
norm. AH tolerance, trust constants, component cap, preconditioner, HVP formula,
coordinates, retraction, and outer energy/gradient thresholds are unchanged.
The old diagonal norm remains available for preconditioning/audit comparisons,
not as the CIAH physical norm.

`test_coupled_structure` compares the matrix-free norm against a dense Riesz
reference and checks the zero covector; `test_ciah` also passes. Molecular
audit 1248698 compares iterative norms with dense diagnostic LLT on actual
keyframes: maximum relative discrepancy is 2.85e-12.

## Short validation results

All molecular runs below stop at an intentional diagnostic iteration cap,
not convergence. The old executable is not rerun.

Four-macro direct comparison:

| Norm | Job | HVPs | Energy, Eh | Physical orbital gradient |
|---|---:|---:|---:|---:|
| Diagonal, corrected interval | 1248684 | 75 | -2181.6191158569 | 0.00825465 |
| Full, corrected interval | 1248694 | 77 | -2181.6195712634 | 0.00982539 |

Energy improves by 0.4554 mEh, but the gradient is larger. This alone is not
evidence of better overall convergence.

Eight-macro comparison using already archived baseline data:

| Implementation | HVPs through step 8 | Step-8 energy, Eh | Physical orbital gradient |
|---|---:|---:|---:|
| Archived original CIAH, 1248081 | 144 | -2181.6244203767 | 0.00982275 |
| Full norm plus corrected keyframe workflow, 1248698 | 149 | -2181.6244495179 | 0.00527183 |

The new eight-step result is approximately 29.14 microhartree lower and has
46.3% smaller physical gradient, for 3.5% more HVPs. Compared with PySCF's
-2181.6244647029 Eh reference, the remaining energy gap decreases from about
44.33 to 15.18 microhartree. The projected infinity gradient is 6.92e-4, but
the last energy change is about 4.95e-5 Eh, so convergence is not reached.
The eight-step historical comparison includes the preceding workflow repairs;
it does not isolate the effect of the norm alone.

1248698 process wall time is 10.27 s and peak RSS 483.8 MiB, including the
dense metric and finite-difference diagnostics. Do not claim a production speedup
or compare this time with a historical full-convergence timing. An audit-free
timing and completed small-system convergence remain necessary before any claim
of total-cost improvement. No iron-porphine long run has been submitted.
