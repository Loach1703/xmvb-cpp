# C12 phase-cost diagnosis, not an optimizer-policy change

## Scope and reproducibility

Only the **current** version was executed. Existing before-change results from
job 1264541 were read, not rerun. The operator fixture is singlet CAS(12,12),
132 covalent structures; the real end-to-end case is LOFLEA, HAO, exact cc-pVDZ,
ISCF=9, Davidson, original input orbitals and numerical settings. This is not
a separate molecular C12 input or a new C14 calculation.

Both jobs below ran on Hanhai22 `anode40`, 32 physical cores, CPU affinity
`60-91`, OpenMP `close/cores`, BLAS threads=1:

- **1264546:** byte-identical c97101a4 production binaries from 1264541. Native
  kernel probes plus the existing phase instrumentation; the SCF completed
  successfully. The subsequent VTune collection failed, making the **batch**
  status FAILED despite a successful application run.
- **1264564:** the same numerical implementation with six finer opt-in phase
  timers and scalar-only CIAH magnitude diagnostics. Release build and **3/3
  relevant regressions passed**; batch status COMPLETED. The SCF completed;
  VTune again rejected collection and the launcher retained that limitation.

No physical formula, convergence criterion, cache budget, tile size, rank
threshold or thread-private numerical workspace was changed. The added phases
reuse the existing `XMVB_PHASE_TIMING`; scalar norms are evaluated only with
the existing CIAH event trace enabled. The build/source hashes and raw logs
are retained in `run-*/` and `evidence.tar.gz`.

Runs 1264562/1264563 were build-maintenance failures, not numerical failures:
the cancellation of 1264562 raced with its start and left zero-byte object
files; 1264563 then failed linking those objects. Only the empty compiler
artifacts were quarantined, and 1264564 rebuilt and tested successfully.
The frozen production binaries remained untouched.

## Does the earlier H/S regression reproduce?

| Native operator timing / s | Saved before, 1264541 | Saved current, 1264541 | Same current binary, 1264546 |
|---|---:|---:|---:|
| H/S cold | 0.149600 | 0.206468 | 0.180126 |
| H/S warm | 0.105439 | 0.127215 | 0.068266 |
| Structure response warm | 0.194061 | 0.189514 | 0.186008 |
| Orbital-direction jet | 0.266877 | 0.285919 | 0.287402 |

The warm H/S value changes substantially with an identical current executable.
Therefore the previous single short probe did **not** establish a stable H/S
regression. It also does not prove a speedup. No old executable was rerun to
obtain these conclusions. The current orbital direction remains about 0.28 s;
the preceding change did not materially reduce its dominant contractions.

The diagnostic-only build's native probes were 0.170260/0.098688 s for cold/warm
H/S and 0.283572 s for the orbital direction. These are a separate instrumented
binary, not additional executions of the unchanged c97101a4 binary.

## End-to-end cost

| LOFLEA | Saved current native, 1264541 | Existing phases, 1264546 | Fine phases, 1264564 |
|---|---:|---:|---:|
| Steps | 4 | 4 | 4 |
| Orbital HVPs | 36 | 36 | 36 |
| Initial energy / Eh | -398.533053417516 | -398.533053417516 | -398.533053417516 |
| Final energy / Eh | -399.106646267331 | -399.106646267330 | -399.106646267331 |
| SCF wall / s | 41.505578 | 41.843810 | 42.298017 |
| Application wall / s | 44.536831 | 44.796305 | 45.591432 |
| GNU-time peak RSS / MiB | 2723.375 | 2676.871 | 2699.379 |

All terminate with `ciah_projected_gradient_tolerance`. The phase runs include
diagnostic overhead and are for attribution, not claimed optimized native
speedups. No added product-square or per-thread large allocation is present.

Within the **36 coupled actions**, the scalar events partition the work:

| Component, fine-phase run | Wall / s |
|---|---:|
| Physical orbital HVP / forcing | 10.315727 |
| CI response and active-integral adjoint | 6.629301 |
| Active adjoint back to orbital coordinates | 0.061417 |
| CI-transport coordinate curvature | 7.578211 |
| Complete coupled-action scope | 24.609468 |

The small difference is metric/coordinate assembly and bookkeeping. Outside
these actions, 11 keyframes cost 6.995613 s and the four macro trial scopes
4.181312 s. The four macro scopes total 37.293619 s. Initial/accepted-point work
outside those scopes and input preparation account for the rest.

Several **nested** measurements locate the costs; do not sum this table:

| Nested component | Calls | Wall / s |
|---|---:|---:|
| Tiled sigma/adjoint response | 88 | 15.055258 |
| Exterior generator pullback | 2488 | 8.164278 |
| Ordinary H/S sigma | 77 | 3.890306 |
| Input panel generation | 812 | 2.525084 |
| Generated CSR patterns | 812 | 1.295417 |
| Product projection | 1412 | 1.797862 |
| Left exterior transforms | 8965 | 1.656182 |
| Supported-transform CSR validation | 1792 | 0.065685 |
| Packing panels | 956 | 0.044726 |
| Small orbital moment pullback | 2488 | 0.028130 |
| Metric solves | 9 | 0.072551 |

The generator timer surrounds the actual function; the moment pullback is
nested inside it. Code inspection shows the generator's spin-string/shear
loops are serial, with `MatrixXd::row().dot(...)` on column-major factor
panels. The roughly 8.14 s excluding its tiny orbital-matrix pullback is not
caused by the small matrix inverse or by CSR checking. These are direct scope
wall measurements, **not** profiler sample percentages.

The fine build additionally instruments the H/S projection callback under
`ci_product_projection`. Its 1412 calls therefore cover more sites than the
earlier 952-call counter; this is expanded measurement coverage, not newly
introduced projection work.

The initial phase run's 165 sigma/response sweeps comprised 3109 batches and
gathered 1,964,120,314 logical doubles, including repeated/zero halo entries.
Aggregate internal stage times were 1.6644 s prepare, 4.9690 s read, 9.3064 s
compute and 2.3953 s write. They are nested within calling phases. These counts
do not mean unique amplitudes or whole panels were regenerated that many times.

## Why removing redundant work did not make the calculation much faster

1. The multi-RHS Gram rearrangement targets the old C14 workspace collapse.
   The C12 selected metric is already cheap; its complete metric solves cost
   only about 0.07 s. This was not the end-to-end bottleneck.
2. Skipping an unsupported target's own diagonal scatter does not eliminate
   its neighbor scatter, dense channel multiplication or reverse contraction.
   In the C12 directional sweep, **304590 targets**, **25 batches** and
   **16872268 gathered doubles** remain. Its fine-phase tile stage costs were
   0.010864 s prepare, 0.061772 s read, 0.096085 s compute, 0.017707 s write.
   Changing only a small diagonal loop cannot remove these dominant operations.
3. CSR reuse removes repeated generation within a sweep, but does not remove
   the serial generator pullbacks in finite CI transport. Patterns are still
   generated 812 times across successive sweeps; packing/checking themselves
   are demonstrably small. Increasing tile capacity is not the inferred cure.

## A mathematically specific next target: stationary-CI transport curvature

The accepted-point constructor in
`src/vbscf/optimization/coupled/ci_transport.cpp` computes, in the equal-weight
convention used here,

$$
g_{c,0}=2w\left(HC_0-MC_0\Lambda_0\right),\qquad
\Lambda_0=C_0^THC_0,\qquad C_0^TMC_0=I,
$$

$$
\lambda_0=M^{-1}g_{c,0}.
$$

Here $H$ and $M$ are the Hamiltonian and overlap operators in the **selected
VB structure space**, $C_0$ contains its relaxed state coefficients, and $w$
is their common state weight. No complete-CAS or closure assumption is made.

Every returned component of `CiTransport::curvature(p,q)` is linear in
$\lambda_0$. This follows directly from `cross_overlap_transport`: the seed
weights are `expand_support(lambda)`, and subsequent operations are linear in
those seeds. Thus at an exactly relaxed selected-space CI point,

$$
HC_0=MC_0\Lambda_0
\quad\Longrightarrow\quad
g_{c,0}=\lambda_0=0
\quad\Longrightarrow\quad
\mathcal K_{\rm transport}(p,q)=0.
$$

This is a **CI-residual-weighted coordinate correction**, not the physical
orbital--CI coupling or CI Hessian. The physical coupled actions need not
vanish when this correction does. Finite-step/keyframe CI-transport gradients
also generally do not vanish and must not be removed on this argument.

Measured over all 36 current coupled directions:

| Transport output norm | Maximum |
|---|---:|
| Orbital from orbital | 1.7191e-11 |
| Orbital from structure | 2.9233e-10 |
| Structure from orbital | 2.1149e-10 |

The maximum ratio of the sum of the two orbital correction norms to the
uncorrected orbital image norm was 1.2729e-10; the corresponding structure
ratio was 7.3846e-10. These are coordinate norms of the observed directions,
not a bound on the whole operator or a replacement for an accuracy contract.

The next mathematical task is to formulate accepted-point stationarity
consistently with the inexact eigenpair contract, eliminating unnecessary
coordinate work without a new cutoff or dropping real CI response. This
round does **not** implement that elimination or claim that the measured
smallness alone licenses it. If generator contractions remain needed, their
shared-panel memory order and independent-output parallelism are the next
implementation target, before any larger cache.

## Profiler availability

`perf` was not installed on the compute-node PATH. The VTune software capture
was rejected because `kernel.yama.ptrace_scope` prevents collection. No
function-sampling or allocation profile is claimed. A candidate Intel
`os-perf` path turned out to be a directory, not an executable; its failed
probe is retained, and the launcher now checks for a regular file. No kernel
settings were changed and no optional profiler was installed.

## Regenerate the machine-readable summary

```bash
python3 tests/mcscf/summarize_c12_cost.py \
  results/diagnostics/c12-cost-20261017/run-1264546 \
  results/diagnostics/c12-cost-20261017/run-1264564 \
  --output results/diagnostics/c12-cost-20261017/summary.json
```

Nested inclusive timers are explicitly labelled in that summary. The previous
C14 convergence result remains untouched; no old version was executed.
