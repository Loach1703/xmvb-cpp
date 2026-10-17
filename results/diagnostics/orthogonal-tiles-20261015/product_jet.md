# Bounded accepted, orbital-direction and reverse product contractions

Parent: `94292e4a`. All C++ builds and numerical runs used Hanhai22 Slurm
compute nodes. No historical executable was rerun. Optimizer strategy,
variational space, integral precision and convergence tolerances are unchanged.
This round completes a memory-ownership boundary; it does **not** establish
a universal runtime improvement over the saved excellent baselines.

## Implementation and exact reuse

- Logical products use signed64-bit addresses, including multiplication before
  packing. The CAS18/18 index test decodes2363904399 without allocating its
 2363904400-entry product space. String and numerical-tile indices remain int.
- Broad-support accepted gradients, orbital directions and backward consumers
  share `product_jet.cpp`. The unbounded accepted Y/residual and directional
  Y/sigma fields and their old dense consumers are deleted.
- Compact accepted Y/residual panels share one64 MiB payload cap across roots,
  inherited from one quarter of the existing256 MiB structure workspace.
  Resident panels supply right readers and residual contractions directly;
  coefficient-only responses do not repeat the accepted exterior transform
  or sigma. Larger states stream primal/tangent images together instead.
  This is a residency decision, not admission to a physical representation.
- The same emitted image supplies selected-space projection and upper
  Cholesky-generator pullback. Sigma and integral adjoint share their channel
  sweep. An accepted identical bra/ket reuses its scatter, with its weight
  applied outside, instead of reading/scattering the same state twice.
  Orbital forward/reverse consumers reuse the prepared integral tangent
  and accumulated small cotangents. No large array is copied per worker.
- Generator beta neighbors borrow compact column views until the next LRU
  access, instead of repeatedly copying their zero-padded halos. Alpha
  reductions use the existing OpenMP team, not a serialized nested team.

Section9 of `article/source_column_hamiltonian.tex` records the accepted/tangent
identities and memory bound. Its20-page PDF was rebuilt and inspected.
Original touched support can itself be complete, and complete-space
inverse/preconditioner consumers are separate remaining owners. A cache-payload
bound is not a whole-process RSS bound.

## Validation

Jobs1262395 and1262404 pass all11 focused tests. Cases include arbitrary
nonclosed selected spaces, both exchange parities, unequal alpha/beta sectors,
multiple states, zero weights, tiny responses and forced128-byte residency
with partial-image eviction. Source and generated-product contraction orders
are compared without relaxing existing tolerances.

Job1262396 runs deterministic constructed exact-Hamiltonian CAS(12,12)
probes with924 strings/spin and1,3,20,132 selected covalent structures.
These are fixed-orbital action/derivative tests, **not molecular SCF runs**.
Rayleigh energies equal their saved printed values; response checksums differ
by at most1.78e-15. HVP central-gradient-difference errors range from5.84e-11
to1.57e-10. Broad orbital-direction reverse consumers take0.00053--0.00055 s
because their forward sweep already supplies the small cotangents; these
numbers are not standalone complete HVP times.

Eight bound OpenMP cores, BLAS1, anode26; native GNU-time/RSS measurements:

| Structures | Fixed-point energy / Eh | Saved best response / s | New response / s | Saved response RSS / MiB | New RSS / MiB |
|---|---:|---:|---:|---:|---:|
|1|-5.6580780473068595|0.07000|0.19074|87.45|48.88|
|3|-5.658608873593056|0.11768|0.23433|86.60|63.29|
|20|-5.658252717096269|0.13730|0.26982|86.88|68.19|
|132|-5.658417026636264|0.13171|0.26441|86.98|70.18|

The single-structure source order is unchanged and its more recent saved
record had48.78 MiB RSS/0.11512 s;48.88 MiB is not a new memory improvement
over that record. The broad cases reduce RSS19--27% against the excellent
saved response runs but remain roughly twice as slow. Do not infer a global
speedup from the lower RSS. Historical and current jobs share the same node
model, but they are separate single runs with other node workloads; this
does not establish a cause for their timing differences.

The supplemental14-orbital/132-structure response reaches298.89 MiB RSS,
versus411.96 MiB in the retained-link-layout record1261842 and438.14 MiB
in the faster record1261729. Its5.57837 s response remains slower than
3.67168/3.33289 s respectively. The checksum agrees within3e-16.

An initial migration had unnecessary dense neighbor copies and nested
OpenMP alpha contraction. The saved current-code experiments1262270,
1262297 and1262346 guided their removal; they were not historical baselines.
At132 structures their response times were0.37062,0.35282 and0.32701 s.
The final0.26441 s is a28.7% repair of that first migration, not an improvement
over the excellent0.13171 s historical response. Keep these comparisons separate.

## Converged molecular regression

Job1262407: exact integrals,32 cores, anode45, unchanged NEO/Davidson and
analytic diagonal preconditioner. All five240 subspaces converge. Steps/HVP
counts are unchanged and final energies differ from saved values by at most
2e-12 Eh. The CSV includes initial energy, application/process time and RSS.

|240 space|Structures|Steps/HVP|Final energy / Eh|SCF / s|RSS / MiB|
|---|---:|---:|---:|---:|---:|
|COV|14|3/20|-343.239211501261|2.12549|1553.34|
|ION(0:1)|294|3/27|-343.407736636064|3.33954|1589.61|
|ION(0:2)|1134|4/32|-343.443415326852|4.71644|1571.35|
|ION(0:3)|1694|4/32|-343.446256564847|4.95238|1591.64|
|FULL|1764|4/30|-343.446302545201|3.78350|1555.87|

The equal-thread/equal-energy excellent FULL reference1259844 was3.39735 s,
1641.07 MiB on anode62. The immediate saved parent1262045 was3.08273 s,
1570.26 MiB on anode21. Current FULL is11.4% slower than the former and22.7%
slower than the latter; its RSS is5.2% lower than the former and0.9% lower
than the latter. All use4 steps/30 HVP, with endpoint differences below2e-12 Eh.
The final allocation is verified as CPUs32--63 on one physical socket.
A prior current-code allocation1262397 split CPUs80--95/176--191 and measured
5.12968 s; it is retained separately, not substituted for the one-socket result
or attributed purely to software. The benchmark now requests one socket and
checks its actual allocation before running numerical work.
The older1.478615 s historical target has a different trajectory/endpoint
(about2e-7 Eh); it remains a speed target, not an equal-solution control.
Do not compare these32-core times with the saved8-core9.94 s run as a
software speedup. Native and opt-in phase runs are retained separately.

## Remaining work

The low-memory accepted/directional ownership loop is implemented, but its
runtime gate is not passed. Next inspect/fuse channel and generator consumers
against saved phase costs; do not enlarge tile budgets, restore unbounded
carriers, change the optimizer, or tune precision to disguise the remaining
cost. Off-diagonal multistate gauge densities still use exact bounded source
contractions. Complete-space inverse/preconditioner workspaces and large
original incidences are not declared solved by this migration.

Machine-readable results: `product_jet12.csv`, `product_jet240.csv`.
`product_jet_evidence.tar.gz` retains input/binary/source hashes, raw logs,
configuration, affinity, GNU-time records and Slurm accounting, excluding
large executable copies retained in the cluster snapshots. A final test-only
unit-oracle update compares coefficient responses with ordinary complete-carrier
sigma/adjoint formulas, not the migrated product jet itself. That test-only
update and final Doxygen cleanup change no timed numerical workload.
