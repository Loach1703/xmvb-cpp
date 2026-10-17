# Residual-balanced CIAH transport: C12 validation

## What changed mathematically

At exact selected-space CI stationarity, the old residual-weighted coordinate
curvature is zero. A numerical CI residual is not exactly zero, so the code
does **not** set it to zero, compare its norm to a cutoff, or conditionally
disable a response. Instead it defines the internal CIAH energy

$$
\widehat E(p,q)=\mathcal E(p,\tau(p,q))
-\langle g_{c,0},\tau(p,q)-(C_0+Qq)\rangle_F.
$$

The actual origin residual $g_{c,0}$ is retained. Because the finite transport
$\tau$ and the affine coefficient map have the same first jet, the
counterterm leaves the origin energy and gradient unchanged and **exactly**
cancels the transport-only second jet, even for a large finite CI residual.
At exact CI stationarity this internal model equals the physical transported
energy for every finite step. For inexact CI it is explicitly a different
internal model, not an assertion that the physical finite-step energy changed.

Both finite keyframe energy and gradient use this same model. Its pullback
uses $g_C-g_{c,0}$ in the existing one metric solve and adds the affine
$Q^Tg_{c,0}$ term back to the CI gradient. The initial CI right-hand side is
not erased. Physical macro acceptance, CI relaxation, energy reporting,
gradients and outer convergence tests remain unchanged.

Full derivation: [CI stationarity and residual-balanced transport](../../../article/ci_stationary_transport.md).

## Implementation boundary

- Keep the physical orbital HVP, orbital-to-CI forcing, CI Hessian/metric
  action, CI-to-orbital adjoint, finite selected-space transport and coarse
  preconditioner.
- Remove the origin residual inverse/adjoint preparation and the separate
  transport-curvature operator, public DTOs, wrapper and fused tile sweep.
  No compatibility/fallback implementation or runtime switch remains.
- Keep the generic selected overlap derivatives and finite transport VJP.
  Selected structures need not form a closed ambient space.
- Add only a scalar energy counterterm and an existing-size CI RHS
  subtraction; no product-space arrays, larger tile budget, or per-thread
  numerical replicas are introduced.

## Verification

Hanhai22 job **1264591**, `anode40`, 32 physical cores, affinity `60-91`, GCC
14.1.0, Release C++23 with the existing optimization flags, BLAS threads=1.
Batch completed with exit 0. All ten edited production/test source files
match the compute-node manifest. Unrelated local RI edits were not transferred.

The initial build job 1264574 stopped on a missed test call to the removed
transport-curvature API, before SCF execution. That test was migrated to the
retained general overlap contractions; no numerical tolerance was changed.

**9/9 regression tests passed** in 3.07 s. Selected-transport checks now cover
both exact stationary roots and deliberately large nonstationary residuals,
single-state and equal-weight two-state cases, on a nonclosed three-structure
subspace. They test:

1. Finite coefficients against independent explicit least squares.
2. Model energy against an independent counterterm definition.
3. Model gradient against energy finite differences at finite steps.
4. Origin HVP against independent affine-coordinate gradient differences.
5. Retention of the nonzero initial CI residual in the nonstationary case.

The physical kernels retain independent numerical checks:

| Independent reference error | Norm |
|---|---:|
| Selected H action | 1.9794e-14 |
| Selected M action | 1.9049e-15 |
| Structure-response image | 5.2861e-16 |
| Structure-response integral adjoint | 1.4979e-15 |
| Orbital-direction HVP finite difference | 5.9108e-11 |

## Native end-to-end comparison

Same original LOFLEA input/guess, CAS(12,12), COV 132, HAO, exact cc-pVDZ,
ISCF=9, Davidson, unchanged tolerances and ITMAX=2000. "C12" denotes the
existing CAS(12,12) fixture and this matching molecular case, not a separate
user-supplied `C12.xmi`.

The saved native baseline is job **1264541 after**, c97101a4. The subsequent
d5df5da3 changes only added opt-in diagnostics. That version's phase run is
used below for attribution, not substituted for a native timing control.
**No old executable or old source was rerun.**

| Metric | Saved native | New native |
|---|---:|---:|
| Initial energy / Eh | -398.533053417516 | -398.533053417516 |
| Final energy / Eh | -399.106646267331 | -399.106646267333 |
| Accepted steps | 4 | 4 |
| Orbital HVPs | 36 | 36 |
| Physical structure responses | 36 | 36 |
| Recorded keyframes | 15 | 15 |
| SCF wall / s | 41.505578 | 33.136628 |
| Application end-to-end / s | 44.536831 | 36.085948 |
| Process peak RSS / MiB (GNU time) | 2723.375 | 2692.195 |
| Mean accounted CPU cores | 26.78 | 27.67 |

Both terminate with `ciah_projected_gradient_tolerance`. The final energy
difference is about 2e-12 Eh. End-to-end wall time is **18.98% lower**, or
**1.234x speedup**. Peak RSS decreases by 31.18 MiB, about 1.14%; this is not
claimed as a new memory-scaling result. GNU-time process peaks and Slurm's
sampled per-task MaxRSS are distinct and are not mixed.

## Where the time went

Separate fine-phase runs, old **1264564** (saved, not rerun) and new
**1264591 loflea-phases**:

| Inclusive diagnostic scope | Saved | New |
|---|---:|---:|
| 36 complete coupled actions / s | 24.609468 | 16.467049 |
| Structure-side scope / s | 14.270589 | 6.593151 |
| Transport-curvature scope / s | 7.574691 | deleted |
| Exterior generator pullbacks | 2488 calls / 8.164278 s | 440 calls / 2.402945 s |
| Physical tiled sigma/adjoint response | 88 calls / 15.055258 s | 88 calls / 14.789951 s |
| Ordinary H/S sigma | 77 calls / 3.890306 s | 77 calls / 3.980793 s |
| Finite transport gradient | 11 calls / 2.709254 s | 11 calls / 2.774114 s |

These scopes overlap and must not be summed. The removed work is the
coordinate-curvature calculation, not fewer physical HVP/response actions.
Finite keyframe pullbacks remain. The largest observed model-versus-physical
keyframe energy difference was 5.6843e-14 Eh across 11 evaluated keyframes;
this measured smallness is **not** an admission rule or the cancellation proof.

The new phase-enabled application took 37.211671 s versus native 36.085948 s;
diagnostic overhead is not included in the stated native speedup.

## Artifacts

- `comparison.csv`: programmatically parsed native energy/counter/time/RSS data.
- `phase-summary.json`: existing saved and new phase data, parsed by
  `tests/mcscf/summarize_c12_cost.py`.
- `evidence.tar.gz`: raw new logs, input, hashes, source manifest, build flags,
  unit details and trace. Executables remain frozen on the cluster under
  `c12-regression-1264591/bin/`.
- Reusable current-version launcher: `tests/mcscf/c12_regression.sbatch`.

No new profiler installation or blocked sampling retry was needed: the
previous investigation had already identified the expensive phase, and this
run verifies its removal with direct scoped timing. No new C14 run was made.
