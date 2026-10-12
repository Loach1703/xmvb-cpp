# FeCl2: tail microsteps and full-metric cost

## Scope

Hanhai22, anode45, one process, 32 bound OpenMP cores, serial BLAS,
original input orbitals, exact integrals, analytic Hessian diagonal,
Davidson, outer gradient tolerance 1e-3 and energy tolerance 1e-7 Eh.
No threshold, microstep budget, HVP, coordinate or keyframe-policy changes.
The existing historical executable was not rerun. No porphine run was submitted.

## Observed tail behavior

Archived job 1248707 used 24 macros and 545 HVPs. Its final eight macros all
ended at `micro_budget`, not only the final six. None of their accepted
microsteps was component-clipped. At macro 24 the predicted dual gradient
fell to 3.938088e-4, then rose to 6.983182e-4 at the budget exit; the source
norm was 6.519534e-4. The last Ritz residual was 6.983181e-4.

Thus a small component limit does not explain this tail. The budget is the
exit reason, but removing it does not address the observed lack of residual
reduction. These logs do not establish a unique cause for ineffective
orbital direction expansion or prove a preconditioner defect.

The inspected installed PySCF sources are `pyscf/soscf/ciah.py:davidson_cc`
and `pyscf/mcscf/newton_casscf.py:update_orb_ci`. PySCF retains `xs`/`ax`,
extends the projected matrix only in its new row/column, and updates the
gradient row at every yield. An accepted approximate microstep changes the
RHS before the AH problem is fully solved. Consequently the sequence is
not convergence of one fixed-RHS eigenproblem. Our implementation similarly
retains images but reconstructs all projected products each time. The
coupled `newton_casscf` implementation is not PySCF's default `mc1step` AH;
these references must not be conflated.

## Accounting and the concrete redundant operation

Job 1248714 added scalar timers and PCG counters without changing numerical
rules. It made 1176 structure-dual-norm calls and 53917 metric actions,
including 52807 CG iterations. The maximum solve took 51 iterations.
These solves consumed 7.065054 s. All AH projection work, including projected
metric certification, consumed 0.312256 s; the generalized eigensolve and
image reconstruction consumed 0.037198 s. The Ritz dual-norm substage consumed
3.438098 s, **nested within** the norm work; do not add it to 7.065054 s.

Every structure metric action previously called the combined H/M action,
discarding H. In direct-CI notation, let D embed the chosen structure vectors
into spin-string products and T denote the exterior orbital transform. Then

$$
M V=D^T T^T T D V,
\qquad
H V=D^T T^T H_{\mathrm{orth}}T D V.
$$

Only forward and adjoint transforms are needed for M. The new `apply_overlap`
entry point removes the Hamiltonian sigma and its carrier entirely on the
direct-CI path. It shares accepted-point data and uses the existing bounded
block width, with no per-thread copies. Pair-based actions retain their
existing exact combined contraction. This change assumes neither complete
structure space nor an inverse of a projected overlap.

## Native results

| Run | Scope | Steps | HVP | Energy / Eh | Wall / s | Peak RSS / MiB |
|---|---|---:|---:|---:|---:|---:|
| 1248707 | archived, before scalar accounting | 24 | 545 | -2181.6244646059 | 24.58 | 498.0 |
| 1248714 | scalar accounting, combined H/M metric action | 22 | 494 | -2181.624464290001 | 22.52 | 540.7 |
| 1248721 | scalar accounting, overlap-only direct-CI metric action | 24 | 544 | -2181.624464442866 | 17.49 | 538.0 |

All three reached the existing stopping conditions. Because trajectories and
HVP counts differ, wall-time ratios are **not** isolated fixed-work kernel
speedups. No numerical-policy change was made, but identical trajectories
were not obtained and their difference is not diagnosed here.

In 1248721, 1290 structure norm solves made 59090 metric actions, taking
1.987587 s. Mean measured metric-action-associated cost fell from 131.04 to
33.64 microseconds, about 74% lower (3.90x). This aggregate includes the CG
vector work and solve overhead; it is not an isolated sigma benchmark.
Maximum certified relative residuals were respectively 1.489628e-8 and
1.489628e-8, below the unchanged square-root machine-epsilon criterion.
AH projection plus spectrum/image work remained 0.357147 s in the new run.

The final eight macros of 1248721 still ended at `micro_budget`, with no
component clipping. This optimization removes a demonstrated cost but does
not solve the excessive HVP count or justify resuming porphine.

## Verification and limitations

`test_ciah`, `test_coupled_structure`, and `test_orthogonal_direct_ci` passed
on the compute node. The latter compares overlap-only to combined action on
three arbitrary sparse structures in a non-invariant selected subspace,
and also checks the incomplete-carrier pair route. Existing dense-reference
tests remain in place. Maximum action difference allowed is 2e-12.

Wall and RSS are native GNU time measurements, not profiler overhead runs.
VTune software sampling in 1248714 could not start because of the compute
node's ptrace restriction. There is no valid hotspot capture, and no sampled
CPU percentage or parallel hotspot attribution is claimed. Stage timers and
action counters establish this explicit redundancy, not a complete profile.

`summarize_ritz_cost.py` parses the scalar traces and TSV automatically.
The next investigation must quantify retained-space Galerkin residuals,
orbital/structure residual components, and the effectiveness of newly admitted
directions across gradient refreshes. Do not increase budgets or loosen norms
to disguise the observed tail.
