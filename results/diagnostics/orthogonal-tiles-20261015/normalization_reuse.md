# Exact reuse of the single-root normalization adjoint

Parent checkpoint: `61e4eb44`. All numerical runs used Hanhai22, anode26,
8 bound OpenMP cores, BLAS 1 and the existing Release build. Historical
executables were not rerun. No precision, optimizer policy or workspace
budget was changed.

## Identity and implementation

Retain the small selected-state gauge frame already computed by the generalized
eigenresponse. For one selected root its coefficient response is exactly
`Q_g = a C`, including in arbitrary nonclosed selected structure spaces.
At fixed accepted active eigenvalue, every quadratic coefficient-dependent
integral or residual-generator adjoint consequently changes by `2a` times its
accepted value. The full-carrier Ritz residual need not vanish.

The exact orbital HVP reuses the retained accepted **active-integral** adjoints
for this contribution. It does not scale the total orbital gradient, which
also contains inactive-only terms. Directional normalization multipliers and
the derivative of the orbital-map pullback remain separate. Coefficient and
multiplier directions now have independent interfaces; no dummy coefficient
array or eager `C A_g` expansion is needed for this case.

The multistate matrix frame is retained, not replaced by the scalar shortcut.
Its general response still includes the necessary cross-state contractions.
The derivation and complexity are recorded in subsection 8.4 of
`article/source_column_hamiltonian.tex`; the rebuilt PDF was visually checked.

## Production gradient/HVP probe, not a molecular SCF timing

Exact constructed 16-AO Hamiltonian, CAS(16,16), one adjacent-pair covalent
structure, 12870 strings per spin. Native phase timers and GNU time are used;
no unavailable sampling-profiler evidence is inferred.

|Quantity|Saved1261411|New1261462|
|---|---:|---:|
|Selected direction / s|55.443474|55.458258|
|Active backward / s|44.843842|20.799125|
|Complete orbital/coupling HVP / s|100.288942|76.258676|
|Extra HVP residual-generator sweeps|2|0|
|Accepted gradient / s|25.885058|25.474356|
|Entire probe including difference gradients / s|177.612053|153.721730|
|Native process peak RSS / MiB|155.3906|155.1133|
|HVP norm|0.583378380827|0.583378380827|
|Central-gradient-difference error|1.7278e-9|1.7278e-9|
|Accepted raw / orthogonal / residual product doubles|0 / 0 / 0|0 / 0 / 0|

The HVP costs 24.0% less than the immediate parent, a 1.315x speedup.
Relative to the saved pre-parity-order HVP of 155.8011 s, the combined changes
give 2.04x; this is not a complete-SCF speedup claim. No additional numerical
array is introduced. The retained source workspace cap remains 64 MiB.
The schedule-derived residual-generator counts exclude the accepted and two
difference-gradient sweeps; raw timer records are retained in the evidence.
Independent allocations can coexist on the node, so short native times remain
exploratory single-run measurements, not exclusive-node scaling results.

## Correctness and complete optimizations

Job1261483 passes 11/11 focused tests. Tests cover the retained gauge-frame
identity with a nonzero Ritz residual, source/product adjoint parity,
independent coefficient/multiplier directions, unequal state weights,
off-diagonal multipliers and integral-direction/no-integral-direction cases.
Production scalar/batch coupled-HVP comparisons and fixed-CI gradient
differences cover arbitrary selected and complete spaces, one/two roots,
with/without an inactive orbital. Single-root tests assert reuse occurs;
two-root tests assert the general frame is used instead.

Job1261467 passes 11/11 tests and nine complete 240 optimizations:

|240 space|Steps|HVP|Final energy / Eh|SCF / s|RSS / MiB|
|---|---:|---:|---:|---:|---:|
|COV|3|20|-343.239211501261|6.561379|1519.16|
|ION(0:1)|3|27|-343.407736636063|8.779741|1518.26|
|ION(0:2)|4|32|-343.443415326853|10.849778|1518.52|
|ION(0:3)|4|32|-343.446256564844|10.971752|1513.98|
|FULL|4|30|-343.446302545201|9.473460|1516.93|

Steps/HVP counts are unchanged; energies agree with the saved immediate-parent
data within 3e-12 Eh. Arbitrary 1/3-structure dense/Davidson runs all take
3 steps and 21 HVP, with consistent final energies. Small-case universal
speedup is not asserted. CSV files retain initial energy and total time too.

Job1261472 passes 8 geometry/coordinate tests and a complete unequal-spin
FeCl2 control, reaching the correct low-energy branch:

|FeCl2 quantity|Saved parent1261419|New1261472|
|---|---:|---:|
|Initial energy / Eh|-2181.617631912267|-2181.617631912267|
|Final energy / Eh|-2181.624464866200|-2181.624464864478|
|Accepted steps / HVP|16 / 271|13 / 235|
|SCF / s|12.086271|10.241470|
|RSS / MiB|416.79|421.70|

Energy agrees within 1.8e-9 Eh. The older correct-metric checkpoint1261264
took 12 steps, 211 HVP and 9.464221 s; current performance does not yet beat
that result. The changed trajectory is recorded, not attributed to rounding
without diagnosis. No convergence-policy change was made here.

## Remaining architectural work

Broad-support product-order accepted gradients and directional consumers still
retain full carrier arrays. Migrate these to the existing generated product-tile
sigma/adjoint graph, preserving oriented generator terms and residual support
outside the selected VB space. This single-root reuse does not finish that
migration or remove the selected-direction cost of roughly 55 s in this probe.

`normalization_reuse_evidence.tar.gz` contains raw logs, hashes, native resource
records, configuration and final focused-test output. Parsed tables are
`normalization_hvp16.csv`, `normalization_240.csv`,
`normalization_selected_240.csv`, and `normalization_fecl2.csv`.
