# Arbitrary selected-space spin metric: diagnosis and first optimization

2026-10-17. Hanhai22/anode42, 32 physical cores, exact Libcint/cc-pVDZ,
LOFLEA HAO CAS(12,12), ISCF=9/CIAH and structure-space Davidson. No optimizer
threshold, orbital/HVP formula, physical structure space or tile memory budget
is changed. No historical executable was rerun.

## What was missing

Arbitrary selected structures already had an exact orthogonal-CI H/S action.
That does not mean the solver enjoyed FULL's complete-sector inverse and
orbital-map cancellations. Those identities require completeness. For large
nonclosed selections, the CI and finite-transport metric preconditioners still
used scalar diagonals and did not remove the Rumer spin-coupling geometry.

The orbital-independent spin Gram `G=W^T W` is occupation-block diagonal for
**any independent subset** of VB structures with definite spatial occupations.
This factorization does not require completeness. The present change extracts
FULL's existing block factor into `structures/overlap/spin_metric.*`, without
duplicating it, and uses a positive block-scaled inverse in selected-space
Davidson and CI transport. The factors are shared through the accepted context.
The nonorthogonal physical `M` is not replaced by `G`, and transport still
solves its original physical equation with the original residual certificate.
FULL retains its exact representation and inner solver.

The block scales use the existing positive Jacobi diagonal estimate, not a new
exact physical-diagonal calculation. One scale per occupation preserves its
spin geometry. ION(0:1) uses **2,002,176 bytes** of factors in **133 blocks**,
maximum block size **132**, not a 5676-by-5676 dense physical Gram. No factors
or integral arrays are copied per thread. Existing 64-MiB factor admission
only controls preconditioner residency, never admission to orthogonal CI.

Derivation and condition bound: `article/source_column_hamiltonian.tex`,
Section 10.7, "Spin-block preconditioning without completeness". The updated
27-page document was compiled with Tectonic and the new section inspected.

## Native complete-run evidence

Control is the already launched **1265093** array, frozen executable
`ec71ee07f7e3677f4fed774cde6130c4d8af392544478275309deacec55bd057`.
Current is array **1265114**, executable
`c51aa5ef1afbbd0ec0259fef2f7f84875a0224bc1a62214626b371238f8e387a`,
frozen by **1265112** after qualification **1265110**: **11/11 tests passed**
in 3.42 s. Subsequent C++ edits only clarify a Doxygen comment about the
Jacobi diagonal estimate; the numerical implementation is the qualified one.
The initial extraction also passed 11/11 in job 1265109 before the final
block-common scaling was qualified. Neither qualification ran old binaries.

The five input files are byte-identical between native arrays. Original guess,
ITMAX=2000, tolerances (energy 1e-7 Eh, projected gradient 1e-3), OpenMP
cores/close placement, BLAS=1 and 64-GiB Slurm allocation are unchanged. Jobs
were submitted as parallel arrays; same-node placement does not imply the
same CPU set or exclusive node occupancy.

| Space | Structures | Control steps/HVP | Current steps/HVP | Control total / s | Current total / s | Current RSS / GiB | Current final energy / Eh |
|---|---:|---:|---:|---:|---:|---:|---:|
| ION(0:1) | 5676 | 4/54 | 4/54 | 416.606497 | 46.759111 | 2.669182 | -399.323427907220 |
| ION(0:2) | 47256 | 6/86 | **5 accepted / 86 total HVP; failed** | 1631.498160 | 173.676211 | 3.066448 | -399.397232509251 |
| FULL | 226512 | 4/40 | 4/40 | 97.036986 | 102.727418 | 3.572872 | -399.417495587055 |
| mixed3 | 3 | 4/38 | 4/38 | 14.640989 | 14.215038 | 2.596825 | -399.050820856823 |
| mixed20 | 20 | 4/38 | 4/38 | 20.185947 | 20.001830 | 2.625374 | -399.095231633343 |

ION(0:1): **8.91x** measured end-to-end acceleration at unchanged outer steps
and HVP count; energy changes by -2.96e-12 Eh. Its control RSS is 2.726921 GiB.
Final CI block actions decrease from 62 to 13. This is lower CI/transport
cost, not a claim of fewer orbital HVPs. FULL is **5.9% slower** in this
single native measurement, with the same steps/actions and energy within
1e-12 Eh; no FULL speedup is claimed. Mixed-subset energies agree within
1.6e-11 Eh. Different spaces are not energy-parity controls for each other.

Total is the application's end-to-end timer, consistently with the existing
C12/C14 report. Full process times (including trailing output/teardown) are
in `comparison.csv`: ION(0:1) 417.15 -> 47.17 s, FULL 107.31 -> 112.47 s.
RSS is GNU-time process peak, not the allocated Slurm memory. The original
ION(0:2) task completed while this report was prepared: 6 steps / 86 HVP,
1631.498160 s, final energy -399.397232509251 Eh. Its energy agrees with
the failed current endpoint to the printed precision, but the current run
is not a formally converged speedup result.

## Attribution, not production timing

One-macrostep ION(0:1) captures **1265104** and **1265126** use GNU gprofng
2.42 plus opt-in scalar phase timers on the respective frozen executables.
Both are deliberately ITMAX=1 and exit 2 at that limit; they are not complete
convergence controls. Their inputs otherwise match. Both collection reports
contain "No errors / No warnings". Native times above exclude profiling.

| First-macro quantity | Before | After |
|---|---:|---:|
| Initial CI solve wall / s | 60.569145 | 2.468535 |
| Relaxed endpoint CI solve wall / s | 51.957516 | 2.104892 |
| Transport solve 1 CG iterations | 389 | 34 |
| Transport solve 2 CG iterations | 1949 | 63 |
| Transport solve 3 CG iterations | 430 | 35 |
| Transport solve 4 CG iterations | 2137 | 70 |

All four new transport true relative residuals remain below sqrt(machine
epsilon), as before. These counts describe each implementation's first-macro
trajectory, not identical RHS vectors. The old capture contains 10 orbital
HVPs costing 2.452979 s, whereas its two relaxed CI solves cost 112.526661 s.
Nested phase scopes are not additive. Sampling resolves sigma, overlap,
exterior transforms and OpenMP waits; sample CPU shares are not stage wall
fractions. Raw `.er` directories remain on Hanhai22 in the corresponding
`c12-ionic-profile-JOBID` directories; scalar logs and profiler reports are
included in this evidence archive.

## Remaining failure: do not count ION(0:2) as converged

The current ION(0:2) run exits 2 with `ciah_macro_energy_rejected` at iteration
6, after five accepted steps. Accepted HVPs total 83; the rejected macro adds
3, for **86** total. Its last accepted energy change was -5.9729748e-6 Eh, so
the normal two-condition convergence test had not yet passed, despite a
projected gradient of 8.49863792e-7.

The rejected trial has source energy -399.39723250925095 and last trial
energy -399.39723250925084 Eh, a computed increase of about 1.1e-13 Eh.
The origin slope is -5.36907243e-13 and curvature 5.36907243e-13.
The unchanged macro policy in `optimization/backends/ciah_step.cpp` requires
strictly positive energy reduction before completing the gradient and letting
the ordinary outer convergence check inspect a candidate. It rejects eight
trials and reports failure. This is a distinct terminal acceptance/convergence
ordering issue, not a measured incorrect H/S action. No threshold was relaxed
and this run is explicitly retained as failed, not relabeled converged.
Resolving that ordering must retain both existing energy and true-gradient
conditions; it must not introduce a gradient-only shortcut.
The original control accepted a final energy decrease of 1.36424205e-12 Eh
and passed its sixth-step convergence check. Thus both endpoints have the
same reported energy, but opposite tiny computed energy changes lead to
different termination statuses under the strict acceptance gate.

## Reproduction from saved evidence

Extract the previous `c12-spaces-20261017/evidence.tar.gz` into a scratch
directory, then extract this directory's `evidence.tar.gz` into the same
directory. The latter overlays the completed control ION logs while retaining
the existing FULL/mixed control records.
Use `tests/mcscf/c12_spaces.py compare SCRATCH/run-1265093 SCRATCH/run-1265114`.
The parser verifies identical inputs, streams scalar output instead of
retaining large structure-weight tables, counts terminal rejected-macro HVPs,
and emits speedup/energy differences only for two completed converged runs.
`comparison.csv` was generated by this parser. The native live directories
are under `/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015` on Hanhai22.
