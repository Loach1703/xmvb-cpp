# Consistent orbital metric and accepted-state producer

## Changes

The metric action now applies the physical tangent map followed by its strict
adjoint, using the same cached accepted-point factors as its squared norm.
The previous expanded action relied on projector identities and separately
computed inverse/transposed factors. Its independent implementation and two
unused AO-by-inactive caches are deleted. There is no new reduced metric,
coordinate whitening, inverse solve, symmetry tolerance or convergence criterion.
The derivation is in `article/orbital_metric_pullback.tex` and its inspected
two-page compiled PDF.

In the complete active-span branch, matrix products are explicitly grouped
to apply the AO operators to narrow tangents/cotangents first. This avoids
forming L*Q_a or Q_a^T*L^T at every action. The old norm's left-associated
L*Q_a expression and the old expanded action's S*Q_a expression both carried
avoidable cubic AO work. The final implementation needs quadratic AO work
times the active/occupied rank instead. The existing accepted-point AO
Cholesky factor is unchanged; no coordinate whitening is restored.

The direct-CI accepted-point producer also stops creating unused raw product
states. NEO's response cost uses the unchanged topology and state dimensions,
not numerical arrays that happen to be resident. This preserves its admission
rule while removing one complete product matrix per selected state.

Remaining dense response producers and broad-support orthogonal/residual
caches are not claimed removed. The source HVP and unrestricted large-space
optimizer still require producer/consumer migration.

## Captured failure and identical-point comparison

All builds/numerical tests ran in Hanhai22 Slurm compute allocations, anode26,
8 bound cores, GCC14.1.0/CMake3.30.5, existing Release O3 flags and BLAS1.
No historical binary was rerun. Job1261206 used the then-current implementation,
no scalar event logging and a diagnostic limit of25 accepted steps. It failed
after13 accepted steps and exported the exact raw orbital table. The AH metric
skew was1.985938952331492e-12, above its unchanged1.874027993755732e-12 bound.

The point-only diagnostic rebuilds orbital geometry without a CI or energy
evaluation. It compares the action with an independent differentiated-feature
Gram matrix. Its dense matrices are diagnostic outputs, never production
optimizer storage. The checkpoint and input hashes are retained.

| Quantity | Saved current-code audit1261210 | Final repaired audit1261264 |
|---|---:|---:|
| Packed / reduced dimension |755/682|755/682|
| Complete action norm |13906.484175711|13906.484175710|
| Complete action skew norm |1.9013342e-10|8.5459663e-12|
| Action minus feature Gram norm |1.1389458e-9|2.5611258e-11|
| Relative skew |1.36723e-14|6.14531e-16|
| Relative feature-Gram gap |8.19003e-14|1.84168e-15|

Skew and feature-Gram gap improve by22.2x and44.5x. The inactive inverse itself
and its discrepancy with the Cholesky inverse factors are unchanged. Thus the
improvement is from reversing the actual feature graph, not changing its
norm, orbitals, conditioning threshold or certificate. A random16-dimensional
projection is also recorded; it is not the original failing AH basis and is
not presented as such.

The saved baseline audit process took4.73 s versus4.91 s for the first repaired
audit1261218. These include the molecular loader, full diagnostic matrices
and independent reference features, not only metric actions. They do not
establish a standalone metric speedup or a large production slowdown.

## Bounded original-input capture, not convergence

| Run | Accepted steps | Recorded HVP | Initial energy Eh | Final energy Eh | SCF / process s | Native RSS MiB | Stop |
|---|---:|---:|---:|---:|---:|---:|---|
|1261206|13|136|-384.380227811199|-395.883611336482|107.605/111.498|2559.18|metric symmetry error|
|1261218|25|265|-384.380227811199|-396.072213699800|200.724/204.748|2567.78|diagnostic iteration limit|

Both application exits are2; the surrounding diagnostic batches completed
after saving the result. The repaired run passes the old failure iteration
and reaches the requested diagnostic limit without that error. It has not
converged and its trajectory differs, so the endpoint energies, HVP counts
and wall times are not a matched convergence speedup. The native failure's
recorded HVP count excludes a possible partially attempted final expansion.
No production iteration limit was changed. The earlier unbounded scalar
diagnostic1261192 was cancelled and is not qualification evidence.

## Complete molecular regressions

Job1261230 passed8/8 focused CTests, including the new cancelling-representative
ray/span adjoint test and a producer test requiring empty raw product state
storage. Its subsequent shared source suite passed9/9 CTests.

The five complete 240 calculations retain3/3/4/4/4 iterations and
20/27/32/32/30 HVPs. Maximum energy difference from saved job1260422 is1.03e-12 Eh.

| Space | Steps | HVP | Saved/current SCF s | Saved/current RSS MiB | Current energy Eh |
|---|---:|---:|---:|---:|---:|---:|
|COV|3|20|6.676/6.568|1528.6/1519.9|-343.239211501262|
|ION(0:1)|3|27|8.761/8.588|1517.8/1524.9|-343.407736636065|
|ION(0:2)|4|32|10.834/10.910|1518.7/1517.7|-343.443415326852|
|ION(0:3)|4|32|10.983/11.160|1515.6/1519.6|-343.446256564846|
|FULL|4|30|9.623/10.062|1525.7/1530.2|-343.446302545202|

These small-space times/RSS are similar, not a universal performance gain.
FULL is about4.6% slower than that saved sample (about5.6% versus saved1170).
No repeats or historical binaries were used to select a favorable time.

The one/three-structure 240 subsets converge in3 steps/21 HVPs for both dense
and Davidson. Paired final energies differ by at most1.03e-12 Eh; SCF times
are6.168/6.359 s for one structure and6.901/6.968 s for three structures.

### Full-AO FeCl2 control

After regrouping the complete-active-span contractions, final job1261264
again passed8/8 CTests and converged the original FeCl2 input with NEO,
Davidson and the analytic Hessian diagonal. It used8 cores:

| Steps | HVP | Initial energy Eh | Final energy Eh | SCF / process s | Native RSS MiB |
|---:|---:|---:|---:|---:|---:|
|12|211|-2181.617631912267|-2181.624464866683|9.464/10.017|407.18|

The final energy differs by5.13e-9 Eh from the saved correct reference
in `current-suite-1249844/summary.csv` (-2181.624464871811 Eh).
That reference used a different node/thread setting; its wall time is not
used as a matched speedup. In particular the much faster old high-energy
branch around-2181.61764 Eh is not used as the correctness baseline.
The final regrouping affects only complete-span contractions; the HAO
point/regressions and the accepted-gradient probe do not execute that branch.
GNU time records758% process CPU, about7.58 effective cores from the8 allocated.

## CAS(16,16) complete accepted-gradient producer

Job1261239 runs `test_keyframe --accepted-state16` through the actual orbital
gradient evaluator, structure-space Davidson, active-integral adjoint and
accepted-point producer. This is a constructed16-AO exact Hamiltonian with
one covalent VB structure, not a molecular SCF or a Hessian/optimizer benchmark.
It has12870 strings per spin. The process reports:

- Complete accepted-gradient evaluation:25.607 s; native peak RSS117.484 MiB.
- Native process CPU:781%, about7.81 effective cores from8 allocated.
- Accepted raw, orthogonal and residual product doubles:all zero.
- The removed raw-product allocation alone:1,325,095,200 bytes (1.234 GiB).

The last value is an exact allocation-size comparison, not a measured old
RSS or speedup. No old allocating executable was rerun. The production-chain
test asserts that process RSS is below that single removed allocation's size;
this is a resource test, never a physical-space admission or approximation.
CAS16 HVP and full convergence remain unqualified.

## Reproduction and next step

Parsed data are `metric_points.csv`, `metric_captures.csv`, `metric_240.csv`,
`metric_selected_240.csv` and `metric_accepted16.csv`. Use
`tests/mcscf/summarize_metric.py`, `summarize_factor_gradient.py` and
`selected_scf.py`; inputs, hashes, raw checkpoints/logs and native resource
records are retained with their job IDs. Batch success is kept distinct from
application exit and convergence.
The full-AO control is additionally parsed in `metric_fecl2.csv`.
Raw evidence is packaged in `metric_consistency_evidence.tar.gz`; unpack it
beside this report to reproduce the relative artifact-directory layout.
Intermediate duplicate point matrices and derived Molden/orbital text files
are omitted from the archive; the original point matrices, final metric,
raw checkpoints, inputs, logs, traces, hashes and resource records are retained.

Next remove the full raw response bundle at its direct-CI producers, including
finite-update gauge coefficient responses. Preserve those responses as
structure coefficient columns through to source-factor adjoints. Then qualify
the complete source HVP producer and its peak RSS at CAS16. Broad product-order
accepted/response arrays and restricted-metric metadata remain later gates.
