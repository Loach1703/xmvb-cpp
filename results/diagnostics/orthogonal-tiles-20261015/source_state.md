# Source-state accepted gradients and orbital/structure directions

## Result and scope

For the source contraction order, accepted states no longer retain a complete
orthogonal amplitude or residual array. The accepted gradient, orbital H/M
direction and complete active-integral HVP pullback use the same exact spin
factors. Directional sigma/coefficients are replaced by a small weighted
generator-gradient direction. The H-image sweep now supplies this adjoint
before its factors are released, avoiding another Hamiltonian traversal.

This completes the source-order derivative consumer chain, not all of VBSCF.
Broad product-order accepted states, their derivative readers, remaining
selected-state producers and restricted-metric metadata still need migration.
No whole-molecule SCF speedup or complete large-active-space optimizer
qualification is inferred from the synthetic active-integral fixtures.

## Execution and qualification

All builds/numerics ran in Hanhai22 Slurm allocations, anode26,8 cores with
explicit srun binding, OMP_PLACES=cores, OMP_PROC_BIND=close and BLAS1.
Final benchmark job1261134 recorded allowed CPUs168-175. GCC14.1.0,
CMake3.30.5 and the existing Release O3 flags were retained. Native clocks,
GNU time and explicit work/storage counters are the measurement sources.
No profiler sample percentage is used as a timing fraction.

Final job1261134 passed9/9 CTests and the numerical, finite-difference and
molecular regressions below. Job1261140 passed9/9 again after a counter-only
cleanup; numerical formulas and benchmark workloads were unchanged.

The first unit job1261081 exposed an error in the new *test oracle*: its sigma
used fixture-combination ordering while its amplitudes used the production
topology's spin-string ordering. The oracle was changed to use the actual
reuse-table ordering; job1261087 then passed9/9. Unqualified overlapping-build
jobs1261103/1261123 were cancelled, their results discarded, and the shared
headers touched to force rebuilding before the isolated qualified runs.
Their statuses are retained; no partial output is treated as evidence.

## Production entry integration

The final audit found that only synthetic tests called
`build_orthogonal_spin_topology()`. All three real active-space forward entries
still supplied only the input's observed one-spin strings. They now complete
the computational one-spin tables and remap the original expansion incidences,
without adding VB structures or independent coefficients. No numerical pair
cache or product amplitude is allocated by the completion itself. The vacuum
topology retains one empty string; its cache-size estimate has zero pair channels.

The new real-input regression uses explicit one/three-structure subsets of
240 and LOFLEA. Job1261164 first failed the new vacuum test because the old
size estimator rejected zero orbitals; this was fixed. Job1261165 passed9/9
unit tests and the five 240 subspaces, but its selected fixture was invalid:
the output's bond notation `27-30` was accidentally used in `$STR`, where it
expands an integer range. The corrected fixture explicitly lists bond endpoints.
Job1261170 then converged for both one-structure solvers and three-structure
Davidson, but exposed a real dense-CI producer omission for three structures:
NEO's CI frame alignment accessed absent selected-state H/M images.

The dense accepted-point producer now supplies those images from its resident
H/M matrices and selected coefficient columns. It does not recompute sigma,
form images for every unselected eigenvector, or change any optimizer criterion.
The selected regressions exercise this contract rather than masking an empty
optional at its consumer. Qualified replacement results are recorded below.

Job1261174 passed9/9 CTests and all four actual 240 selected-space SCFs:

| Structures | CI solver | Steps | HVP | Final energy Eh | SCF / total s | Native RSS MiB |
|---|---|---:|---:|---:|---:|---:|
|1|Davidson|3|21|-342.941087950298|6.159/8.542|1517.13|
|1|dense|3|21|-342.941087950298|6.295/8.634|1526.40|
|3|Davidson|3|21|-343.040331009666|6.762/9.099|1513.73|
|3|dense|3|21|-343.040331009664|7.012/9.323|1521.04|

The paired initial energies agree exactly at printed precision, and final
energies differ by at most2e-12 Eh. These are complete molecular SCFs, not
kernel times. No saved pre-change run exists for these exact selected decks,
so no before/after whole-SCF speedup is assigned to them. Their purpose is
production routing and the dense/Davidson accepted-state contract.

The completion does not by itself remove the raw accepted-product state
still allocated by `finalize_active_space_second_order_context()`, or the
corresponding response producer. Those are explicit next milestones; an
empty source-order amplitude/residual cache must not be confused with a
fully product-array-free production optimizer.

### Open failure: real CAS(12,12) selected LOFLEA

Job1261175 passed9/9 unit tests, then its one-covalent-structure Davidson SCF
accepted24 steps (274 completed recorded orbital HVPs), reducing the energy
from-384.380227811200 to-396.139660977277 Eh. It then failed with
`CIAH metric action is not symmetric`. SCF time was199.890 s, process wall
203.863 s and native peak RSS2547.60 MiB. The remaining three runs were not
started because the batch stopped on failure. None is reported as converged.
The top-level exception report prints0 iterations; the CSV instead retains
the24 actual accepted iterations from the saved trajectory. A possible
partially attempted next HVP is not included in the recorded-action total.

This case has no nonzero CI tangent: the offending AH metric action is the
orbital chart/physical metric, not a selected-structure CI metric solve. That
localizes the failing subsystem but does not prove why symmetry was lost.
The source contraction's derivative comparisons do not certify this operator
at this evolved molecular point. The unchanged symmetry threshold is retained.
Job1261192 enabled scalar diagnostic events and final raw-orbital export, with
the same numerical operators and criteria. It passed the original failure
iteration on a different trajectory (at iteration24 the actual reduction was
0.003324599268467 instead of0.003150693078623 Eh), then continued beyond49
accepted steps (50 recorded) without that failure. It was cancelled at7m12s rather than
treated as a successful reproduction or extended blindly. It did not reach
the final exporter, so no exact failing-point raw orbital table was obtained.
The new exception text can report existing skew/tolerance quantities on
failure but those numbers have not yet been captured for the native failure.
No routine actions or new stopping criterion were added. Capture the precise
native failing point before choosing a metric repair; the full CAS12 molecular
qualification remains open.

## Representation and implementation

- `source_direction.cpp` generates delta-H Y, H delta-Y and the left-generator
  H Y projection directly from factor panels. Both sides of the changing
  orbital representation are retained; the selected VB map is not assumed closed.
- `source_image.cpp` contracts raw H-image factors into weighted right-spin
  panels. Common exchange parity is reconstructed on the left by its exact
  half-weighted transpose, not inferred from equal alpha/beta dimensions.
  Only after the complete image is consumed are the oriented generator links
  contracted. Additional right panels may replay H sources, but never repeat
  the structure projection.
- Image and adjoint workspaces partition one64 MiB numerical budget. The two
  additional projection-derivative slabs and provider scratch are included.
  Data is shared across workers; only small orbital reductions are thread-local.
  Original support vectors, returned tensors and immutable metadata are separate.
- `source_response.cpp` also reuses Hq and HY for both response adjoints.
  Hq provides H/M structure images; HY needs no structure projection.
  The operation-count planner was updated from12 to6 times U*r*d^2 for this
  fused request, retaining all source-pair moment terms. No fitted threshold
  or physical-space admission was introduced.
- State-average weights, off-diagonal multipliers and the differentiated
  Cholesky/integral pullback are unchanged. Standalone gradient-direction calls
  without an already produced H-image adjoint still compute the required
  residual terms; those terms are not dropped.

The equations and memory contracts are in Sections8.4-8.6 of
`article/source_column_hamiltonian.tex` and the rebuilt PDF. These are actions
and covectors, not explicitly constructed Hessian blocks or high-order RDMs.

## Direction plus HVP timing

The fused direction stage includes its generator-image adjoint, so compare
the **sum** of direction and HVP times, not HVP time alone. The product
reference uses the existing cached-product directional/residual order with
the same integral-adjoint planner; it is also matrix-free, not a dense Hessian.
Reference processes are current-code independent oracles, not historical binaries.

| Active space / structures | Unfused source1261090 s | Final source1261134 s | Current product oracle s | Final native RSS MiB |
|---|---:|---:|---:|---:|
|12/12 /1|0.38598|0.28649|0.28545|55.64|
|14/14 /2|7.72892|4.13100|6.75156|76.81|

CAS14 is46.6% faster than the first uncached source implementation and38.8%
faster than the current product oracle. CAS12 is essentially tied with the
product oracle, not a demonstrated universal speedup.

For CAS12/CAS14, source H visits are128/384 columns, exactly4*r for r=32/96;
the unfused graph's three residual contractions plus images require12*r
analytically. This is an internal channel-work reduction, not fewer optimizer
HVP calls. The cached generator direction holds only144/196 doubles.
Accepted/directional full-product counters are zero in both native runs.

The final source H/M direction and HVP stages separately took0.22881/0.05769 s
and2.98799/1.14301 s. Their accumulated numerical workspace bounds are
33,136,224/55,691,424 bytes, both below64 MiB.

## Matched accepted preparation and structure-response process

Same benchmark input, source coefficients, state energy, node,8 cores and
native command as saved jobs1261034 (CAS12) and1261052 (CAS14). A process
includes setup, accepted preparation and cold/warm response calls, not SCF.
Historical executables were not rerun.

| Active space / structures | Accepted prep saved/current s | Response warm saved/current s | Process wall saved/current s | RSS saved/current MiB |
|---|---:|---:|---:|---:|
|12/12 /1|0.17573/0.09340|0.07000/0.11979|0.36/0.35|87.45/49.25|
|14/14 /2|2.97863/1.71197|2.83700/1.72163|8.83/5.31|293.02/78.59|

The CAS12 isolated response warm time remains slower than its saved sample;
it is not hidden or called an optimization. The whole process is comparable
and uses less memory. A one-structure horizontal CI tangent is zero: this
scalar response case tests algebra/plumbing, not a nontrivial CI search.
CAS14 uses two covalent structures and a nonzero response projected orthogonal
to the accepted vector in the structure metric. It shows a39.9% process-time
reduction and73.2% lower native peak RSS.

## Numerical evidence and CAS16

New unit tests compare source and independent product gradients/directions
for nonclosed selected spaces, both exchange parities and no common parity,
unequal spin counts, unequal state weights, cross-state gauge terms, left/right
orbital motion and the full connection diagnostic. Cached orbital directions
plus coefficient responses are also checked together. Existing finite-gradient
differences and HVP symmetry checks remain enabled.

CAS12/CAS14 accepted-gradient errors are6.36e-15/7.53e-15; H/M-direction errors
are4.34e-19/3.91e-19, and HVP errors9.30e-18/2.03e-17. The independent
unconstrained-generator/residual kernels remain separately tested.

CAS16 accepted-gradient preparation takes22.69 s and its native process peaks
at149.88 MiB, including two shifted-action finite-difference evaluations.
No accepted product array is retained. CAS16 has12870 strings
per spin; its two amplitude/residual arrays alone would need2,650,190,400 bytes
(2.47 GiB), before any HVP temporary. This is an exact removed-allocation
count, not a measured old RSS or timing. No CAS16 product oracle was run.
The Lagrangian directional finite-difference errors are6.49e-11 (CAS12),
3.91e-11 (CAS14), and3.13e-11 (CAS16). CAS16 H/M directions, HVP and complete
orbital convergence remain unmeasured; their CSV fields are blank.

## Complete 240 regression

The production-entry regression in job1261170 uses exact libcint, NEO/Davidson, analytic Hessian diagonal,
projected-gradient threshold1e-3 and energy-change threshold1e-7 Eh. Compare
to saved job1260422 on the same node/eight cores, without rerunning it.

| Space | Steps | HVP | Saved/current SCF s | Saved/current RSS MiB | Current energy Eh |
|---|---:|---:|---:|---:|---:|
|COV|3|20|6.676/6.517|1528.6/1525.1|-343.239211501262|
|ION(0:1)|3|27|8.761/8.863|1517.8/1507.2|-343.407736636063|
|ION(0:2)|4|32|10.834/10.710|1518.7/1518.2|-343.443415326852|
|ION(0:3)|4|32|10.983/10.811|1515.6/1530.2|-343.446256564845|
|FULL|4|30|9.623/9.527|1525.7/1517.9|-343.446302545201|

Steps/HVP are unchanged and the maximum energy difference is2.96e-12 Eh.
Whole SCF time/RSS remains essentially unchanged for these broad-support
cases. There is no claim that this source-consumer milestone removes the
remaining broad-product or restricted-metric allocations.

Parsed data: `source_state.csv`, `source_state_response.csv`,
`source_state_240.csv`, `source_selected_240.csv`, `source_selected_loflea.csv`,
accounting and raw hashes/logs. Reproduce parsing with
`tests/mcscf/summarize_factor_gradient.py`. GNU-time native RSS is used for
runtime memory; the batch/compiler RSS and allocation-only CPU rows are not
application-memory or useful-core measurements.
