# COV CI: metric-aware Davidson corrections and Hamiltonian-only tiles

## Why 132 structures needed 99 actions

The old COV solve used a generalized Davidson problem in the original Rumer
structure coefficients, with the self-spin-product Jacobi approximation.
FULL instead used its complete-sector isometric adapter and congruent carrier
preconditioner. FULL's 27 actions and COV's 99 were therefore not evidence
that smaller CI spaces intrinsically require more work.

The selected COV metric has condition number 295,190 at the original orbital
point. Scaling by the existing approximate overlap diagonal gives 297,497:
diagonal normalization does not remove the correlations between Rumer
columns. At the new trial point these figures are about 385,722 and 388,726.
These are measured metric conditions, not a claim that the overlap diagonal
is the entire old energy-dependent Davidson denominator.

The important mathematical distinction is covector versus tangent. The Ritz
residual r = Hc - E Mc is a covector. Its metric representative is M^-1 r.
With a Ritz basis V, Galerkin orthogonality V^T r = 0 gives
V^T M (M^-1 r) = 0. A raw coordinate diagonal does not preserve this relation
in general. The full selected metric inverse maps the residual into a
correction in the appropriate geometry, without changing the eigenproblem
or enlarging its variational space.

This is a metric preconditioner, not an exact inverse of H - E M. No
energy-dependent parameter or denominator floor was tuned. Davidson still
operates in the **132-dimensional VB structure space**. String-product
arrays remain internal contraction intermediates, not solver variables.

## Implementation and eliminated work

1. Use the bounded selected metric already supported by `PreparedMetric`.
   Factor it once per CI orbital point and solve Mz = r by two triangular
   solves. The correction needs no exterior transform or sigma action.
2. Share that same selected M with the Davidson H/M callback. In the product
   contraction order, compute H alone through generated tiles and supply Mv
   from the small selected matrix. The alpha adjoint, beta sampling and
   structure projection for M are no longer repeated per CI action.
3. Check matrix-storage admission before constructing optional metric
   transforms. This is the existing metric workspace bound, not a condition
   that enables/disables orthogonal CI. Large selected spaces keep their
   bounded generalized action and existing inexpensive preconditioner.

The ordinary source-factor order retains its fused H/M action. Its H and M
consumers are not claimed to have been split in this round. The new
Hamiltonian-only API returns the same H on both orders; it omits the
redundant product-order overlap projection where the CI adapter owns M.

The 132-by-132 M and its factor each occupy 139,392 bytes, about 0.133 MiB.
They are shared per CI problem, not copied per thread or action. Construction
and the small factor remain inside the existing 8 MiB metric allowance;
there is no full string-pair matrix or additional carrier cache. No dense
structure Hamiltonian is formed. The CI callbacks, coordinates, HVP,
orbital preconditioner and stopping tolerances are otherwise unchanged.

Section 10.6 of `article/source_column_hamiltonian.tex` records the covector,
Galerkin and basis-congruence derivations and the storage/contraction scope.
Its 26-page PDF was compiled and the new section visually inspected.

## Native one-step qualification

Same LOFLEA geometry/original HAO guess, exact libcint cc-pVDZ, CAS(12,12),
132 COV structures, ISCF=9 NEO/CIAH, structure Davidson, analytic orbital
Hessian diagonal, gradient 1e-3 and energy 1e-7 Eh. Eight actual OpenMP cores
within one LLC, single-threaded BLAS; 16 CPUs and 48 GiB allocated.

| Version/job | Node | SCF (s) | Process wall (s) | Peak RSS (GiB) | HVP/response | Final CI actions |
|---|---|---:|---:|---:|---:|---:|
| Saved remote 94693457 / 1262622 | anode50 | 53.464025 | 58.06 | 2.612 | 6/6 | 99 |
| Previous current f979999e / 1263242 | anode04 | 30.096030 | 34.55 | 2.531 | 6/6 | 99 |
| Selected metric preconditioner / 1263304 | anode50 | 16.886731 | 21.41 | 2.522 | 6/6 | 20 |
| Metric plus H-only product projection / 1263331 | anode50 | 16.341363 | 20.88 | 2.497 | 6/6 | 20 |

The final one-step process wall is **2.78 times faster (64.0% less time)**
than the saved remote on the same anode50. Remote binaries were never
executed in this round. Relative to the previous current, the observed
reduction is 39.6%, but that comparison crosses nodes. There are no repeated
timing averages. The two new native runs isolate the metric-only and H-only
implementation stages, not repeated samples of one binary. Their small
0.53 s difference alone is not a robust standalone H-only speedup estimate.

All one-step runs have 10 microsteps, three keyframes, six HVPs and six
structure responses. Initial energy is -398.533053417516 Eh in every case.
After one step, the previous current gives -399.074297727156 Eh and the final
version -399.074110938218 Eh: they are different optimization trajectories,
not identical-point timing samples. This 0.187 mEh difference is not hidden
as a final energy accuracy claim. Full convergence was checked below.
Final one-step CI stationarity is 8.73378158e-8 Eh, below the unchanged 1e-7
budget, and relative residual is 3.17766554e-12.

## Full convergence and historical energy check

Current-only job 1263333 ran from the same original input without the one-step
limit. It converged normally; exit status is zero.

| Result | Accepted steps | HVP/response | Final energy (Eh) | Gradient norm | Process wall (s) | RSS (GiB) |
|---|---:|---:|---:|---:|---:|---:|
| Final COV / 1263333 | 4 | 36/36 | -399.106646267332 | 1.16471042e-5 | 71.48 | 2.510 |

SCF time is 67.355929 s, program end-to-end time 71.236984 s. The last energy
change is 1.4211e-10 Eh. The final CI needs nine H/M actions and its residual
is 6.15272696e-8 Eh. No convergence threshold or extra certification was added.

Existing exact-integral converged LOFLEA output 1257927 gives
-399.106646267326 Eh, also four accepted steps and 36 HVPs. Initial energy
is identical; the final difference is about **6e-12 Eh**. Its saved input
differs only by STR keyword case and ISCF=7, overridden by its explicit NEO
command. Geometry, basis, HAO coefficients, active space and tolerances match.
That older run used 32 cores on anode42, so it is used here as a numerical
reference, not a matched performance comparator. It was not rerun.

## Tests and attribution

Final build/test job 1263323 passed 12/12 tests on Hanhai22. New tests compare
Hamiltonian-only with joint H/M actions on covalent and unequal-spin
nonclosed selections. Selected metric tests use independent occupied-minor
references and verify both inverse residuals and the congruence relation
P (P^T M P)^-1 P^T r = M^-1 r. Existing CI transport, root, keyframe,
source/product order and derivative checks continue to pass.

| Explicit phase | Previous current / 1263226 | Metric only / 1263305 | Final / 1263332 |
|---|---:|---:|---:|
| Sigma calls | 205 | 45 | 45 |
| Sigma time (s) | 12.998060 | 2.730940 | 2.630086 |
| Left-transform calls | 11,616 | 4,544 | 3,899 |
| Left-transform time (s) | 5.425954 | 1.619358 | 1.273106 |
| Orbital HVP time (s) | 2.067272 | see CSV | 1.988619 |
| Initial CI solve (s) | 10.067935 | 2.122145 | 2.023248 |
| Trial CI solve (s) | 9.457978 | see raw log | 1.516798 |

These nested phase clocks are attribution, not additive CPU sampling shares.
The previous-current phase run is on anode04; both new phase runs are on
anode50. Native runs, not phase runs, establish the end-to-end figures.
The principal gain is the reduction of CI expansions, with an additional
strict removal of the repeated product-order M projection. Remaining cost is
now distributed among coupled CIAH directions, keyframes and AO/setup work;
the small metric triangular solve is no longer an evidenced bottleneck.

## Evidence and source scope

The comparison, native, converged and phase CSVs are code-generated. The
archive retains native/phase/convergence stdout, traces, GNU-time, Slurm,
input/binary hashes, CPU/affinity, the final build and all test output, plus
the scoped patch and source hashes. The failed first prototype build 1263283
was a missing header, fixed before any workload ran; successful builds
1263297, 1263309 and 1263323 passed all 12 checks. No queued task remains.

Unrelated local RI edits and `.gitignore` changes were preserved, neither
included in this commit nor transferred as part of this work.
