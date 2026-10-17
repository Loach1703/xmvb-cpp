# Exact source-column H/M contraction for selected VB spaces

## Implemented change

Parent revision: `f531a61b`. All runs used Hanhai22 Slurm compute nodes;
historical executables were not rerun.

The orthogonal-CI representation and structure-space Davidson are unchanged.
There are now two exact contraction orders within that representation:

1. Generated product tiles, suited to broad source support.
2. Original source-column factors, suited to few touched source columns.

The source order separates alpha/beta factors throughout the Hamiltonian
action and samples only original touched products before the final structure
contraction. It handles arbitrary nonclosed selected spaces, unequal spin
dimensions and mixed exchange sectors. Certified common exchange parity also
eliminates duplicate original input products; the output keeps both halves.
Selected-row exterior pullbacks retain every causal intermediate edge needed
by the requested output rows.

The operation-count planner is constructed once per immutable action. It
compares channel GEMMs, spin links, exterior shears and final sampling.
Neither memory nor a fitted threshold changes the physical space. Source and
channel panels jointly obey a 64 MiB explicit numerical-workspace bound.
The allocated slabs are shared and reused across all panels and RHS; no
per-thread carrier, global transformed product, or U-squared pair cache is
introduced.

Theory: [source_column_hamiltonian.md](../../../article/source_column_hamiltonian.md).
This scaling change currently covers fixed-point H/M actions, not every
integral adjoint, orbital gradient or HVP consumer.

## Final validation

Final job **1260459**, anode26, 8 bound cores, GCC14.1.0/CMake3.30.5, Release
`-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`, BLAS1. All eight CTests passed.
The new tests exercise both source and channel streaming, partial tiles,
multiple RHS, both parities, even diagonal products, metric-only source actions, asymmetric selected
spaces, unequal spins, and selected-row adjoints against complete adjoints.

CAS(12,12) tests use 1, 3, 20 and 132 covalent structures. The explicit
source order is also tested where the planner chooses product tiles. Against
the independently expanded full-carrier reference, maximum H-image error is
2.54e-14 and maximum M-image error is 3.56e-15. Native Rayleigh energies
differ from the saved baseline by at most 3e-15 Eh.

## Matched native CAS(12,12) H/M kernel

Saved baseline1260276 and final1260459 both used anode26 with 8 cores.
The fixture, width1 and orbital/integral values are identical. These are
fixed-orbital synthetic H/M actions, not molecular SCF runtimes. RSS is the
whole native process peak; the independent reference runs are separate.

| Selected structures | Saved warm ms | New warm ms | Saved RSS MiB | New RSS MiB | Chosen order |
|---:|---:|---:|---:|---:|---|
| 1 | 59.800 | 36.386 | 68.277 | 50.496 | source, 32 columns |
| 3 | 72.845 | 63.262 | 68.266 | 50.715 | source, 64 columns |
| 20 | 90.931 | 87.032 | 70.391 | 70.633 | product tiles |
| 132 | 113.490 | 85.115 | 71.660 | 72.098 | product tiles |

The final source-order gains are **1.64x** and **1.15x** for 1/3
structures, with approximately 26% lower process RSS. Product-order work
for 20/132 is unchanged. Its short timing differences are not assigned to
the source algorithm: in particular the saved baseline132 independent
reference process already had an 88.678 ms warm action, close to the new
85.115 ms. Raw comparisons, including these differences, are retained in
`source_comparison.csv`.

The preceding workspace-reuse verification1260411 measured35.909/49.415 ms
for1/3 structures. The final cleanup removes unused Hamiltonian allocations
from metric-only actions and adds an independent metric check. Its combined
H/M order is unchanged. Both sets are retained; the table uses the final
native measurements, rather than selecting the faster3-structure sample.

Forcing source order for 20/132 structures costs 142/394 ms despite giving
the same H/M images. The exact contraction-order choice avoids that excess
work; no legacy representation or numerical approximation is involved.

## CAS(14,14) single covalent structure

There are **3432 strings per spin**, 128 original touched products and
64 representative source columns. The native action uses 22 source panels
and requires no 3432-by-3432 product array.

- Native cold/warm action: **0.267 / 0.262 s**.
- Action construction, including planning: **0.071 s**.
- Native process peak RSS: **52.152 MiB**.
- Rayleigh energy: **-5.7651827186932545 Eh** for the synthetic fixture.
- Separate full-carrier reference process: H error **2.67e-15**, M error
  **8.89e-16**; reference evaluation **2.141 s**, source warm action in
  that same process **0.257 s**.

The full-carrier reference includes its own integral/exterior preparation;
it is an independent numerical oracle, not the historical production tiled
baseline. No historical14-orbital speedup is claimed.

## Allocation diagnosis and correction

The initial source implementation passed parity but was too slow. Saved
candidate1260396 has the same representative source reduction but repeatedly
allocates/free its large channel panels. Its14-orbital native process had
**723689 minor faults** and **6.55 CPU seconds in system mode**.

Final1260411 allocates each source/channel slab once and uses full-row views
for partial panels. Minor faults fell to **21528**, system CPU to **0.09 s**,
and warm action from about **0.8 s to0.266 s**. These are saved candidate
results; the preceding binary was not rerun. Earlier unquotiented candidates
had approximately1.6 s warm actions.

Final native phase timing for one14-orbital action reports approximately
55 ms scatter,66 ms channel mixing,44 ms gather,67 ms selected pullback and
0.5 ms sampling. These explicit phase clocks are not profiler sample shares.
VTune was absent at the checked compute-node path; that probe is retained.
There is no sampling capture and no profiler percentage claim.

## Molecular convergence regression

Job **1260422**, anode26, **8 cores**, exact libcint, NEO/Davidson with
analytic diagonal preconditioning. All five240 subspaces converged. The
saved comparison1260261 used anode15/**32 cores**; it is used only for
numerical and action-count checks, not a speedup claim.

| Space | Structures | Steps | HVP | Final energy Eh | New SCF s | New RSS MiB |
|---|---:|---:|---:|---:|---:|---:|
| COV | 14 | 3 | 20 | -343.239211501261 | 6.676 | 1528.605 |
| ION(0:1) | 294 | 3 | 27 | -343.407736636066 | 8.761 | 1517.754 |
| ION(0:2) | 1134 | 4 | 32 | -343.443415326851 | 10.834 | 1518.676 |
| ION(0:3) | 1694 | 4 | 32 | -343.446256564847 | 10.983 | 1515.570 |
| FULL | 1764 | 4 | 30 | -343.446302545201 | 9.623 | 1525.727 |

Steps and HVP counts match the saved runs; energy differences are at most
5e-12 Eh. These complete240 subspaces retain product order under the
arithmetic planner, so they validate the surrounding pipeline rather than
demonstrate a source-factor SCF speedup. There were no tolerance or optimizer
strategy changes.
