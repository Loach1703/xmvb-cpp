# Work-planned parity placement in the integral tangent

Parent checkpoint: `dd0ffc3d`. The preceding exact two-family merge did not
improve total HVP time. This round changes the contraction order responsible
for repeated directional-source generation; it does not enlarge the budget.

## Exact identity and implementation

For a certified common exchange projector, a self-adjoint physical integral
variation obeys

$$
\langle P_\eta\delta Z,\delta\mathcal H Z\rangle_F
=\langle P_\eta Z,\delta\mathcal H\delta Z\rangle_F.
$$

The two source-count pairs are `(4r,r)` and `(2r,2r)` after the exact tangent
merge. They have the same source-pair moment count. The first has less leading
channel-Gram work but may repeatedly regenerate a large directional bra. The
second can retain the smaller accepted bra in the existing numerical slabs.

`SpinSources::integral_tangent()` compares both exact forms using
`mcscf::detail::factor_integral_work()`. The estimate includes the actual
planned provider loading, source/channel replay and leading Gram work; it
uses the same slab planner as the numerical kernel. No molecule-dependent
constant, amplitude cutoff, closure assumption, convergence threshold or
physical-space admission test is added. Without a certified common parity,
the original bilinear form is used. Hamiltonians are **not** moved through
oriented residual generators, where this identity would be invalid.

The identity is documented in `article/source_column_hamiltonian.tex`,
Eq. 56. The PDF was rebuilt and visually checked, with no overfull or
undefined-reference warnings.

## Matched production HVP probe

Hanhai22, anode26, 8 bound OpenMP cores, BLAS 1, identical Release flags,
same exact constructed 16-AO Hamiltonian, CAS(16,16), one adjacent-pair covalent
structure, 12870 strings per spin. This is a production-gradient/HVP/difference
probe, **not** a converged molecular SCF benchmark. Both runs use native
opt-in phase timers, not an external profiler. The previous column reuses
saved diagnostic job 1261387; no historical executable is rerun.

|Quantity|Saved1261387|New1261411|
|---|---:|---:|
|Integral-adjoint bra / ket sources|512 / 128|256 / 256|
|Left / right / channel slab widths|14 / 8 / 12|256 / 1 / 20|
|Numerical workspace / bytes|66934752|65626480|
|Scheduled bra column loads|106496|256|
|Provider generation / s|52.9624|0.6491|
|Integral-adjoint timed subtotal / s|75.6640|20.1406|
|Selected structure direction / s|55.2698|55.4435|
|Active backward / s|100.5300|44.8438|
|Complete orbital/coupling HVP / s|155.8011|100.2889|
|Full probe including both difference gradients / s|231.6461|177.6121|
|Full process peak RSS / MiB|155.3164|155.3906|
|HVP norm|0.583378380827|0.583378380827|
|Central-gradient-difference error|1.7278e-9|1.7278e-9|
|Accepted raw / orthogonal / residual product doubles|0 / 0 / 0|0 / 0 / 0|

The HVP is 1.55 times faster; the entire gradient/HVP/difference probe is
1.30 times faster. The workspace cap remains **64 MiB**, and no large array
is replicated per worker. GNU time reports 790% CPU use for the new probe.
RSS is its process maximum, not the size of one slab or a Slurm polling value.

The source-load numbers are schedule-derived counts from the fixed kernel
loops, not independent allocation-profiler observations. For the saved case,
512 columns are replayed across 16 ket panels and 13 moment passes:
`512 * 16 * 13 = 106496`. The new accepted bra is loaded once. Provider
generation time, scatter, moments and Gram time are directly instrumented.
Nested general phase totals are not treated as additive wall-time fractions.

Raw numerical values and kernel records are parsed automatically into
`tangent_order_hvp16.csv` by `tests/mcscf/summarize_metric.py`.

## Numerical and complete-SCF verification

- 1261395: 9/9 focused tests, including both parity signs, mixed anchors,
  unequal spins and full-carrier integral oracles. A mock expensive loader
  forces the alternative exact order, verifying that it is exercised rather
  than merely present in the source.
- 1261413: 9/9 focused tests and 9 complete 240 SCF regressions. Every case has
  unchanged steps and HVP count and agrees in energy with the saved reference
  within 2.0e-12 Eh.

|240 space|Steps|HVP|Final energy / Eh|SCF / s|RSS / MiB|
|---|---:|---:|---:|---:|---:|
|COV|3|20|-343.239211501262|6.576620|1515.90|
|ION(0:1)|3|27|-343.407736636064|8.722837|1519.79|
|ION(0:2)|4|32|-343.443415326851|10.859024|1513.05|
|ION(0:3)|4|32|-343.446256564845|11.022926|1518.61|
|FULL|4|30|-343.446302545200|9.389819|1510.36|

Arbitrary selected 1/3-structure dense/Davidson optimizations all take 3 steps
and 21 HVP, with mutually consistent final energies. Small-system speedups
are not asserted. See `tangent_order_240.csv` and
`tangent_order_selected_240.csv`.

Job 1261419 also passes 8 geometry/coordinate tests and a full FeCl2 control:

|FeCl2 quantity|Saved correct metric checkpoint1261264|Current1261419|
|---|---:|---:|
|Initial energy / Eh|-2181.617631912267|-2181.617631912267|
|Final energy / Eh|-2181.624464866683|-2181.624464866200|
|Accepted steps / HVP|12 / 211|16 / 271|
|SCF / s|9.464221|12.086271|
|RSS / MiB|407.18|416.79|

The energy agrees within 4.9e-10 Eh, but this control is slower and takes more
steps than that older checkpoint. Its saved data predates both the response
producer migration and the two-family merge, so it is not a matched immediate
parent comparison. FeCl2 has unequal alpha/beta spaces and no common exchange
quotient; the new projector-placement order does not apply there. This is
remaining regression evidence, not a universal optimization claim.

## Remaining work

The selected-direction phase still costs 55 s, and two normalization-response
residual adjoints each cost about 12 s. For a known single-state gauge response
`Q_g = a C`, its coefficient-dependent adjoints scale exactly as `2a` times
the accepted adjoints. Reuse must be based on that upstream algebraic frame,
not numerical detection of a nearly parallel vector. The multiplier and
orbital-map derivative terms must remain.

Broad-support accepted and directional product consumers still need streaming
migration. This source-order improvement does not finish that requirement.
No high-order RDM or orbital Hessian is explicitly formed by this change.

Raw logs, hashes, input decks, resource records and build configuration are in
`tangent_order_evidence.tar.gz`. The profiler failures are separately retained
in `tangent_merge_evidence.tar.gz`: gprofng is unavailable and VTune collection
was denied by ptrace policy. No CPU sampling percentages are claimed.
