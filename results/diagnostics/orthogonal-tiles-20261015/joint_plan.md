# Joint source/channel scheduling at a fixed memory budget

## Mathematical cause and implementation

The saved CAS16 residual kernel selected left/right/channel widths1/1/105.
Its128 ket panels each consumed the256-column bra four times:131072 source
columns were regenerated, and moment GEMMs repeatedly rewrote large outputs
with rank-one updates. The failure was not the exact residual-factor formula.

The new shared `mcscf/ci/factor_plan.cpp` searches feasible channel and ket
widths and maximizes the bra width for each pair. It compares nominal scalar
work and element touches from source regeneration, spin-link operations,
channel zero fills and output updates. Leading invariant GEMM arithmetic is
not counted again. Source work comes from writes, diagonal scales, exterior
shear operations and oriented-generator work, not fitted molecule parameters.
There is no new environment option, accuracy threshold or physical admission.

Both integral and generator adjoints consume this plan; the two old
channel-first planners were deleted. Partial-channel residual sweeps currently
use one ket column to remain consistent with the shared forward workspace's
channel extent. This restriction remains an optimization opportunity.

Single-replacement links are now generated directly in packed-channel order
from opposite occupations. Scatter/gather visits only the selected interval
using binary ranges in the existing arrays. No sorting copy, full channel
offset table, duplicate integral data or per-thread spin array is allocated.
Matching-link visits fall from quadratic to linear in channel-panel count;
binary lookup work remains and is included in the schedule estimate.

Theory and nominal-work/resource formulas are in Section8.4 of
`article/source_column_hamiltonian.tex` and the rebuilt PDF. The estimate is
not a guarantee of hardware-optimal wall time; native results qualify it.

## Native CAS16 result

Hanhai22 anode26,8 bound cores, GCC14.1.0/CMake3.30.5, Release
`-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`, BLAS1. Same physical single
covalent structure, orbital/integral fixture and64 MiB workspace budget.
These are kernel times, not SCF or complete orbital-HVP times.

| Version/job | Ordinary generator s | Residual generator s | Native peak RSS MiB |
|---|---:|---:|---:|
|Saved1260873|0.151948|174.370516|118.80|
|Joint plan1260921|0.152790|15.250428|118.51|
|Joint plan + channel ranges1260948|0.142810|12.124799|118.52|

Final residual-kernel gain is **14.38x**. Old candidates were not rerun.
The final plan uses widths256/1/20 and66,983,392 bytes of explicit numerical
workspace, below64 MiB. Phase counters confirm **256 columns in one bra load**,
matching the analytical schedule. RSS includes immutable links, exterior
metadata and runtime allocations beyond that workspace.

The independent native process took25.06 s for cold/warm generator and
residual calls, versus349.58 s saved. CPU accounting stayed near eight cores
(793% versus797%); voluntary context switches fell from1482885 to33998.
These are whole-process counters, not sampled attribution or useful-SIMD
utilization percentages.

## Numerical and larger-scale verification

Numerical/application job1260948 passed9/9 CTests, including arbitrary
nonclosed selected spaces, unequal spins, parity, directional gradients/HVP,
multistate gauge terms and small-workspace streaming. Final supplemental job
1260991 passed all9 CTests and also exercises every individual channel against its independent oriented
one-body generator, and check actual provider loads against planned counts.

For CAS16, a separate phase/reference process checks a nonsymmetric generator
direction without allocating a12870-by-12870 carrier. It computes the
Hamiltonian bilinear from the integral adjoint paired with h/V, and the norm
bilinear from the number-operator trace, then compares with the residual
generator contraction. Scalar errors are6.94e-18 and8.33e-17. This is an
independent contraction-order identity, not a full-matrix reference; unmeasured
full-matrix error fields remain blank in the CSV.

## Saved smaller kernels and complete 240 regression

Same anode26/eight-core conditions. Native warm times are retained without
choosing the faster reference-process samples.

| Kernel | CAS | Saved s | Current s | Saved/current RSS MiB |
|---|---|---:|---:|---:|
|Integral adjoint|12/12|0.05516|0.03379|73.50/47.75|
|Integral adjoint|14/14|0.43506|0.41364|82.49/82.79|
|Directional integral adjoint|12/12|0.06656|0.05556|74.91/51.12|
|Directional integral adjoint|14/14|0.63638|0.60394|75.22/75.23|
|Residual generator|12/12|0.05438|0.05782|53.78/48.93|
|Residual generator|14/14|0.60423|0.53578|61.66/61.62|

The small CAS12 residual warm sample is slower, while its cold sample and
reference-process factor sample are faster; all raw observations are retained.
Single short runs do not establish uniform speedup or a precise noise model.

Job1260948 also ran complete 240 COV/ION/FULL optimization, exact libcint,
NEO/Davidson with the same Hessian-diagonal preconditioner and convergence
criteria (projected gradient1e-3, energy change1e-7 Eh).

| Space | Steps | HVP | Saved/current SCF s | Saved/current RSS MiB | Current energy Eh |
|---|---:|---:|---:|---:|---:|
|COV|3|20|6.676/6.661|1528.6/1507.3|-343.239211501262|
|ION(0:1)|3|27|8.761/8.794|1517.8/1519.8|-343.407736636063|
|ION(0:2)|4|32|10.834/10.925|1518.7/1519.5|-343.443415326853|
|ION(0:3)|4|32|10.983/10.954|1515.6/1539.1|-343.446256564844|
|FULL|4|30|9.623/9.421|1525.7/1517.5|-343.446302545202|

Steps/HVP counts are unchanged and maximum final-energy difference is2.96e-12
Eh. Whole SCF time/RSS is essentially unchanged in these broad-support cases.
The14.38x claim applies only to the measured CAS16 residual kernel.

## Remaining scope

Accepted orthogonal/residual carriers, other generator/HVP consumers,
selected-state producers and restricted metric metadata still need migration.
The new schedule does not itself remove their full-product allocations.
Continue with fused integral/generator moments and complete consumer migration;
do not put the new residual kernel on top of an already-cached residual and
thereby repeat Hamiltonian work.

Data and reproducible parsing: `joint_plan_generators.csv`,
`joint_plan_integrals.csv`, `joint_plan_240.csv`, raw logs/hashes and
`tests/mcscf/summarize_factor_gradient.py`. Measurements are native clocks,
GNU time and explicit phase counters, not unavailable profiler sampling.
