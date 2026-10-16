# Fused, bounded orthogonal-CI adjoints

## Result

The requested fusion, beta parallelism and bounded group-projection accumulation
are implemented without changing the optimizer, physical representation or
precision. Numerical validation passes. The performance gate against excellent
saved baselines is **not yet passed**. In particular, the 240 convergence run
does not show a material end-to-end speed improvement.

Only new code was run. All old results below are retained artifacts, not reruns.
Compilation and all C++ execution were Hanhai22 Slurm compute jobs.

## Changes

- The generated response panel supplies both the transposed accepted-residual
  contraction and the energy-multiplier contraction in one generator traversal.
  Transposition is applied before the upper-triangle mask.
- Alpha/beta contractions share the existing OpenMP team and small orbital-size
  reductions. Beta spectator dots are parallel and SIMD-reduced. Links discarded
  by every requested orientation are excluded before sorting/gathering.
- Immutable accepted column views are borrowed directly; transient LRU views
  are copied into shared, panel-bounded spectator batches before parallel use.
  No worker owns a product array or structure-size numerical replica.
- An accepted reader schedules its whole requested column batch once, rather
  than launching a team separately for every source panel.
- Reverse projection accumulates linear contributions in bounded group tiles.
  Only tile eviction/finalization traverses the existing group-to-structure
  incidence. This removes the per-emitted-panel whole-structure scan and
  repeated immediate structure scatter. It does not assume selected-space
  closure or discard small coefficients.

The existing 64 MiB cache budget and transform extent are unchanged. Projection
tiles share a 64 MiB payload cap across their output consumers; their row extent
is additionally bounded by one existing numerical panel. An output tile may
coincide with a complete small support, but complete-product storage is never
required. Large supports evict exact partial cotangents. This is not a physical
representation/admission branch. Metadata and process allocator retention are
not represented by the numerical payload cap, so it is not a whole-process RSS
bound. This work does not eliminate the original selected-expansion metadata.

## Validation

- Build/CTest job1262521, anode20: 11/11 focused tests pass. Coverage includes
  unequal spins, exchange parities, multistate/zero weights, arbitrary selected
  subspaces, tiny responses, and forced 128-byte residency/eviction.
- Added generator-specific tests compare scaled/transposed fused requests with
  independent complete-carrier and resident-column formulas. Both immutable
  and invalidated borrowed-column ownership are tested inside/outside an
  existing OpenMP team.
- Native action/gradient/HVP probes: job1262522, anode26, eight actual OpenMP
  cores, BLAS1. The allocation reserves16 CPUs to select eight cores sharing
  one LLC; timings do not represent16-thread execution.
- Separate accuracy-only job1262527 runs the **current binary's** complete-carrier
  oracle for H/M and coefficient-response actions at1,3,20,132 structures.
  Oracle memory/times are excluded from native production measurements.

Maximum absolute norm errors across those four cases:

| Quantity | Error |
|---|---:|
| H action vs complete-carrier oracle |1.95e-14|
| M action vs complete-carrier oracle |2.02e-15|
| Response H image vs oracle |8.89e-16|
| Response integral adjoint vs oracle |2.78e-15|
| Gradient vs central energy difference |6.50e-11|
| HVP vs central gradient difference |1.57e-10|

## CAS(12,12) native probes

These are constructed fixed-orbital probes with924 strings/spin, not converged
molecular SCF calculations. Time is the complete warm coefficient response;
RSS is the native process peak measured before invoking a separate oracle.

| Structures | Excellent saved / s | Previous low-memory / s | Current / s | Excellent saved RSS / MiB | Previous RSS / MiB | Current RSS / MiB |
|---|---:|---:|---:|---:|---:|---:|
|1|0.070003|0.190743|0.123111|87.45|48.88|49.21|
|3|0.117684|0.234335|0.141791|86.60|63.29|60.28|
|20|0.137304|0.269820|0.167314|86.88|68.19|63.99|
|132|0.131711|0.264413|0.146996|86.98|70.18|70.83|

The excellent single-source reference1261034 used more accepted storage. The
more recent source-order reference1261793 was0.115122 s/48.78 MiB. This round
does not change that source-order algorithm and does not beat its saved time.
For3--132 structures, current RSS is19--30% below the excellent references,
but response time remains12--22% slower. Fixed-point energies equal the saved
printed values; checksum differences are at most1.78e-15.

The14-orbital/132-structure supplement is3.741987 s/281.04 MiB. The previous
low-memory record was5.578374 s/298.89 MiB; the excellent saved record1261729
was3.332887 s/438.14 MiB. The current response checksum differs by less than
3e-16. Native internal cold-response timing is used for this one-sweep case.

### Placement matters; do not attribute the entire change to the algorithm

The previous12-orbital run1262396 used CPUs167--174, crossing two LLC domains.
The excellent saved1261723 used152--159, one LLC. The current native run uses
168--175, one LLC. In the **same final executable/allocation**, the132-structure
control on167--174 takes0.256320 s, versus0.146996 s on168--175. Thus a large
part of the apparent0.264→0.147 improvement is placement, not software.
`cache-affinity.jsonl` preserves the actual masks and LLC grouping.

Opt-in phase records (separate from native timings) show generator calls falling
from120 to90 across accepted+cold+warm sweeps. Their cumulative time changes
from0.100337 to0.044461 s, but the placement difference also applies. Sigma's
17 batches,304590 targets, and tangent source-element count13816132 are
unchanged. No reduction of that underlying channel work is claimed.

## Converged240 regression

Job1262523, anode45,32 cores on one socket, CPUs32--63. Exact integrals,
NEO/Davidson, analytic diagonal; gradient1e-3 and energy1e-7 Eh, unchanged.
All five spaces converge with the saved step/HVP counts. Energy differences
from the saved table are at most2.05e-12 Eh.

| Space | Structures | Steps/HVP | Final energy / Eh | SCF / s | Native peak RSS / MiB |
|---|---:|---:|---:|---:|---:|
|COV|14|3/20|-343.239211501259|2.127704|1561.47|
|ION(0:1)|294|3/27|-343.407736636064|3.271556|1596.15|
|ION(0:2)|1134|4/32|-343.443415326852|4.690135|1564.80|
|ION(0:3)|1694|4/32|-343.446256564845|5.177818|1584.65|
|FULL|1764|4/30|-343.446302545200|3.734777|1628.86|

| Saved/current FULL record | Steps/HVP | SCF / s | Peak RSS / MiB |
|---|---:|---:|---:|
|Excellent same-energy1259844, anode62|4/30|3.397351|1641.07|
|Faster immediate-parent1262045, anode21|4/30|3.082727|1570.26|
|Previous low-memory1262407, anode45|4/30|3.783504|1555.87|
|Current1262523, anode45|4/30|3.734777|1628.86|

Current FULL is9.9% slower than1259844 and21.2% slower than1262045. The1.3%
change versus1262407 is not evidence of a material overall speedup. Current
native RSS is higher than1262407; do not claim a universal RSS reduction.
The current separate phase run peaks at1572.27 MiB, illustrating that one
native peak cannot identify an allocation owner. Saved timings are single runs
on the same CPU model, not controlled old/new repeats on identical node load.
The older1.478615 s/12-HVP record has a different trajectory/endpoint and is
not substituted for these same-energy4/30 controls.

## Remaining cost

The fused generator and bounded projection are numerically qualified, but the
excellent-baseline performance gate remains open. The remaining12-orbital
channel sweep still does the same scatter/GEMM work, and generated halo reading
and panel/image pullbacks remain important. Full240 does not benefit enough
from this CI-local change; avoid claiming that a faster isolated response
establishes faster molecular convergence. Further work should target measured
remaining traversal/dispatch rather than enlarged budgets or changed precision.
The orbital-direction Hamiltonian and integral-direction Hamiltonian still use
two separate sigma sweeps; this round has not implemented their joint channel
sweep. It must not be counted as an already-delivered fusion.

Machine-readable results: `fused_product12.csv`, `fused_product240.csv`.
`fused_product_evidence.tar.gz` contains native logs, separate parity/phase
runs, CTest output, input/binary/source hashes, affinity, compiler configuration,
and Slurm accounting. Cluster artifacts retain the exact executable copies;
the local archive excludes executable binaries. `fused_product_sources.tar.gz`
contains the modified production/test sources and benchmark launchers.
