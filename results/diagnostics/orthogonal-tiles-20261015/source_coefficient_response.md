# Structure-space response producers

## Change

Parent checkpoint: `fe251298`. Direct-CI directional and pure response-adjoint
interfaces now accept structure coefficient columns instead of expanded raw
spin-product matrices. Exact-HVP relaxed and gauge responses borrow their
existing columns. The non-direct pair consumers keep their separate required
representation; no compatibility overload remains in `StructureAction`.

The exact differentiated representation remains

$$
\delta Y=\mathcal G_\alpha(\kappa)Y+
Y\mathcal G_\beta(\kappa)^T+T_\alpha W\delta C T_\beta^T.
$$

Only the last term uses the original touched support. The two orbital terms
are not projected onto it. Selected-space closure is not assumed. All
normalization and full state-multiplier terms remain, including off-diagonal
multistate terms. Pure adjoints share the Hq/HY factor sweep but omit unused
structure H/M output projections. Product-order adjoints expand one response
column locally rather than retaining an upstream multi-state raw bundle.

Theory: `article/source_column_hamiltonian.tex`, Section 8.2; compiled PDF checked.

## Environment and evidence

Hanhai22, anode26, 8 bound OpenMP cores, BLAS 1, Release
`-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`, GCC 14.1.0, CMake 3.30.5.
Native GNU time supplies RSS;
Slurm step RSS is reported as zero and is not used. Existing old results are
read from saved CSV files; no historical executable was rerun.

- 1261329: build and 9/9 focused tests, including factor/product HVP and
  multistate response parity, unequal spin strings and arbitrary supports.
- 1261339: 9/9 focused tests plus 9 complete SCF regressions.
- 1261340: 9/9 focused tests and the production-gradient/HVP CAS16/16 probe.

## Complete SCF regressions

All 240 calculations use exact integrals, NEO/Davidson, analytic-diagonal
preconditioning and unchanged stopping criteria. The previous column is saved
job 1261230, on the same node with 8 cores. Time is SCF iteration wall time,
excluding initial integral generation and final output.

|240 space|Steps|HVP|Final energy / Eh|Previous / s|New / s|Peak RSS / MiB|
|---|---:|---:|---:|---:|---:|---:|
|COV|3|20|-343.239211501259|6.568444|6.548304|1544.19|
|ION(0:1)|3|27|-343.407736636064|8.588487|8.742364|1516.88|
|ION(0:2)|4|32|-343.443415326852|10.909604|10.879681|1518.69|
|ION(0:3)|4|32|-343.446256564846|11.160497|11.043311|1523.82|
|FULL|4|30|-343.446302545201|10.061604|9.722080|1520.04|

Steps and HVP counts are unchanged. Maximum energy difference from 1261230
is 3.02e-12 Eh. Small-system times are essentially unchanged; no general
speedup is claimed for this producer migration.

Arbitrary selected 1/3-structure spaces also converge with both dense and
Davidson CI: 3 steps, 21 HVP in every case. Final energies agree with the saved
results within 2.1e-12 Eh and with each other within 2.1e-12 Eh. Their full
numbers are in `source_coefficient_selected_240.csv`.

## Production CAS16/16 HVP and memory gate

This is a **constructed exact 16-AO Hamiltonian**, one actual adjacent-pair
covalent VB structure, and 12870 strings per spin. It is not a molecular
convergence or wall-time benchmark. The actual orbital-gradient evaluator,
structure Davidson, accepted-point producer, exact orbital/coupling HVP and
two relaxed-gradient evaluations are invoked, not only low-level factor tests.

|Quantity|Measured value|
|---|---:|
|Accepted raw / orthogonal / residual product doubles|0 / 0 / 0|
|Accepted gradient|25.731710 s|
|One complete orbital/coupling HVP|153.490708 s|
|HVP norm|0.583378380827|
|Central-gradient-difference error, step 1e-4|1.7278e-9|
|Full probe wall time, including both difference points|230.04 s|
|Full probe peak RSS|155.6953 MiB|
|CPU use|792%, approximately 7.92 effective cores|

The test also asserts zero raw selected-response reconstruction time and empty
accepted product caches after HVP. A single complete raw product matrix would
contain 165636900 doubles, or 1325095200 bytes (1.234 GiB). This is an allocation
size, **not** a measured old process RSS. No old allocating executable was run.
The complete gradient/HVP/difference probe stays below that one-array floor.

The 153.49 s HVP still needs substantial work. No same-call production-HVP
baseline exists at this size, so no HVP speedup factor is inferred. Next steps
are stage attribution, exact sharing of coefficient/orbital tangent factors,
and streaming the remaining broad-support accepted and directional consumers.

Parsed data: `source_coefficient_240.csv`,
`source_coefficient_selected_240.csv`, `source_coefficient_hvp16.csv`.
Raw logs, inputs, affinity, hashes, native resource records and build flags
are retained in the three source-cover directories and factor-gradient1261339
and archived in `source_coefficient_response_evidence.tar.gz`.
