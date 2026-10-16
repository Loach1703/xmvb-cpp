# Generated product response and compact symbolic columns

Parent checkpoint: `c0a7a368`. All numerical work and builds ran in Hanhai22
Slurm jobs on anode26, with eight bound OpenMP workers and BLAS threads one.
Existing executables were not rerun. Native timing uses the existing Release
flags. Scalar phase timing is enabled in the paired 14-orbital captures;
these are not sampling-profiler results. Other jobs share this node, so the
small timing differences are exploratory single-run measurements.

## Algebra and ownership

The selected-coefficient response now uses the same generated product graph
as the structure H/M action. With accepted Y and R=(H-E)Y borrowed, the
generator response is

    2 [ D(Hq,Y) - E D(q,Y) + D(q,R)^T ].

The input columns supply the last two terms once; the emitted Hq columns
supply both the first term and structure projection before release. The same
channel scatter supplies the integral adjoint. There is no full response q,
Hq, response-residual, or projection array, and no repeated HY action.
The upper-triangular physical Cholesky cotangent uses the previously proved
support certificate. General oriented generator contractions are not
silently projected, and transpose precedes upper-triangle selection.

Compact column caches retain values on the exact symbolic support, with
shared sorted row indices across RHSs and input/output panels. Structural
zeros and tiny amplitudes inside that support are retained. Sigma gathers
directly into compact rows; projection unpacks one bounded panel. The two
cache limits are unchanged at 64 MiB each for this fixture. Full support uses
implicit indices. There is no admission test on the physical representation.

Exterior cache misses now share the existing sigma OpenMP team through
disjoint-row/column tasks. They do not start serialized nested teams, copy
floating-point carriers per worker, or change any shear's coefficient order.

## Cache-thrashing diagnosis and measured repair

The completed native process in job1261623 records both cold and warm calls;
its subsequent reference process was cancelled. Preserve that distinction:
the native results exist, but this job supplies no completed reference check.
Job1261729 makes one native response call and uses a separate reference
process. Its independent complete-carrier differential validates the new
response, rather than calling the migrated implementation as its own oracle.

Both are the same deterministic CAS(14,14), 132 selected covalent structures,
3432 strings per spin, one coefficient-response direction, exchange parity +1.

| Quantity | Saved 1261623 | New 1261729 |
|---|---:|---:|
| Accepted preparation / s | 3.18775 | 3.17182 |
| Cold complete response / s | 293.93278 | 3.33289 |
| Warm response / s | 285.18316 | not run |
| Native process peak RSS / MiB | 1661.664 | 438.145 |
| Input columns generated as cache panels | 11024 | 54 |
| Input evictions | 10986 | 0 |
| Output evictions | 12183 | 0 |
| Peak input cache / bytes | 66772992 | 29064108 |
| Peak output cache / bytes | 66470976 | 29064108 |
| First sweep input-reader / s | 131.52139 | 0.95669 |
| First sweep channel compute / s | 1.73511 | 1.43049 |
| First sweep output consumer / s | 160.07246 | 0.32939 |
| Sigma batches | 317 | 317 |
| Generated channel targets | 3886025 | 3886025 |

The complete cold response is about 88 times faster. This is elimination
of exact-zero padded residency and repeated exterior/projection work, not
an 88-fold reduction in Hamiltonian channel arithmetic or a claimed full-SCF
speedup. Native RSS is the process high-water mark, not only cache payload;
no allocation-profiler attribution is inferred from its reduction.

New independent errors: structure H/M image norm 1.11e-15 and combined
active-integral/overlap adjoint norm 6.45e-16. The native checksum changes
from -0.077403701794113239 to -0.077403701794114044.

## Smaller cases and tests

Jobs1261723 and1261729 each pass all 11 focused tests. Cases cover arbitrary
nonclosed selected structure spaces, unequal spins, parity 0/+1/-1, weighted
multiple roots, no-image adjoints, and response amplitudes of order 1e-180.
Compact-row sigma is checked against the complete output before projection;
inside-team exterior primal, adjoint, inverse, directional and row-cache-miss
operations are compared with their ordinary-call results.

CAS(12,12) response warm times, compared with saved pre-compaction1261622:

| Structures | Order | Saved / s | New / s | New RSS / MiB |
|---:|---|---:|---:|---:|
| 1 | source factors | 0.11745 | 0.12515 | 49.30 |
| 3 | product tiles | 0.12347 | 0.11768 | 86.60 |
| 20 | product tiles | 0.15993 | 0.13730 | 86.88 |
| 132 | product tiles | 0.15078 | 0.13171 | 86.98 |

The single-structure factor order was not changed by column compaction and
does not show a speed improvement. Its older 0.07000 s record1261034 remains
a separate performance target. Small-case process RSS is mostly unchanged;
at 132 structures each column cache falls to 2721900 bytes while accepted
carrier arrays remain resident.

Machine-readable results are in `product_response_compact.csv`. Raw logs,
source/binary hashes, CPU model, affinity, CMake flags and accounting are
retained with the job artifacts. Complete 240 selected/COV/ION/FULL SCFs and
the FeCl2 control qualify the final integrated revision separately below.

## Complete optimizations

Job1261730 passes all 11 tests and completes nine 240 optimizations with the
new code: selected1/3 with both dense and Davidson, and COV/ION/FULL. No
optimizer strategy or stopping threshold is changed. All five family runs
preserve their step and HVP counts, with saved-baseline energy differences
below 3e-12 Eh. The following immediate saved parent is1261633; historical
data are reused, not rerun.

| 240 space | Steps/HVP | Saved SCF / s | New SCF / s | New RSS / MiB |
|---|---:|---:|---:|---:|
| COV | 3/20 | 6.6673 | 6.3685 | 1518.59 |
| ION(0:1) | 3/27 | 8.9514 | 8.5650 | 1499.03 |
| ION(0:2) | 4/32 | 10.7280 | 10.0897 | 1521.66 |
| ION(0:3) | 4/32 | 11.1247 | 10.5578 | 1524.63 |
| FULL | 4/30 | 9.6221 | 8.9645 | 1514.09 |

Selected1/3 still take3 steps and21 HVP; dense/Davidson final energies differ
by at most2e-12 Eh. Those four SCF times are5.6895,5.7766,6.5747,7.1357 s.
Process peak RSS is1516.02--1524.49 MiB. Original input preparation and
integral storage remain in those complete processes.

Job1261731 passes the additional eight geometry/metric/optimizer tests and
converges FeCl2. This is **not** a universal performance win:

| FeCl2 quantity | Saved1261657 | New1261731 |
|---|---:|---:|
| Accepted steps | 16 | 16 |
| Orbital HVP | 269 | 276 |
| Final energy / Eh | -2181.624464849574 | -2181.624464870515 |
| SCF / s | 11.420771 | 13.207529 |
| Peak RSS / MiB | 406.242 | 415.742 |

The energy difference is2.09e-8 Eh under unchanged outer criteria, but the
action count and wall time regress. Earlier saved correct controls were also
faster (1261545:15 steps/249 HVP,11.1858 s). Do not explain these trajectory
differences as rounding without analysis or infer a per-HVP speedup from the
14-orbital response probe. They remain a separate optimization target.

The automated tables are `product_response_240.csv`,
`product_response_selected_240.csv`, and `product_response_fecl2.csv`.

## Remaining ownership boundary

Broad-support accepted Y/R arrays and orbital-direction product arrays have
not been removed. Original touched-support scratch may itself fill a complete
product grid. This change therefore does not certify the entire VBSCF chain
as free of product-sized allocations. The next migration must address those
accepted and orbital-direction consumers together, not delete their producer
and force downstream recomputation. No selected-space closure or approximate
wavefunction representation is introduced.
