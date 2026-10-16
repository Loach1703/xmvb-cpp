# Single-structure endpoint cost and exact source-jet reductions

## Findings before changing the algorithm

The existing LOFLEA SINGLE input has 185 Cartesian basis functions, 28 inactive
pairs and CAS(12,12), with one six-singlet-bond HAO VB structure. Its original
support touches only 64 spin products. Exact spin-exchange reduction yields
32 source factors; the one-spin computational tables each have 924 strings.
These dimensions must not be confused with 924-by-924 numerical amplitudes:
the accepted and directional source path does not allocate that product array.
A one-dimensional selected CI problem also does not make its correlated
many-electron wave function a single spin product.

Current-at-entry phase capture 1263412 establishes that the initial CI solve
costs only 0.045365 s. The initial forward structure setup costs 0.009808 s.
They cannot explain the 36.51 s complete optimization. More significant costs
are six first-macro HVPs (2.231320 s), including their selected/integral
direction (1.294953 s), and the analytic Hessian diagonal (2.254031 s).
Inside that diagonal, density contraction costs 2.052229 s. The exact
preconditioner uses 52 density channels: 28 inactive plus 24 active source-span
channels, in batches of eight. This is not a large CI eigensolver problem.
The molecule-static exact AO graph has 73,962,904 unique integrals; initial
integral generation/input already costs about 4.3 s in saved native outputs.

Phase clocks are explicit nested scopes, not additive sampling percentages.
VTune capture 1263414 was denied by the compute node's ptrace policy; its
stderr is retained. No function-level profiler percentages are inferred.

## Mathematical removals

Section 10.7 of `article/source_column_hamiltonian.tex` derives both changes.

1. **Residual VJP fusion.** The Hamiltonian source sweep already emits the
   factors needed by its oriented generator adjoint. Emit the identity source
   with coefficient minus the state energy once per source panel. This forms
   the residual adjoint in the same sweep. Delete the additional overlap
   generator contraction and its separate physical-factor construction from
   both the orbital-direction and coefficient-response consumers. Returned H/M
   images are unchanged; the energy-direction multiplier remains separate.
2. **One-column image contraction.** For one selected VB column, obtain its
   scalar H/M directional images by contracting accepted one-electron,
   Coulomb-kernel and generator adjoints with their integral/Cholesky
   directions. Account for the actual scalar CI coefficient and state weight:
   the API returns delta H times C, not its quadratic expectation. No exterior
   output adjoint, output generator or structure projection is necessary for
   those scalars. The full directional integral/generator adjoints, CI
   normalization transport and orbital pullback remain. A requested connection
   diagnostic still computes its separate uncontracted images/norm; it is not
   enabled in native production runs.

This is an exact elimination of repeated operations, not a different variational
space, frozen response, new convergence test, rank truncation, larger tile
budget or modified preconditioner. Multicolumn arbitrary selections retain the
exact restricted-space images and benefit from the residual fusion where their
source order is selected. Broad product-order consumers are unchanged.

## Native endpoint comparison

All use the original identical input/HAO guess, exact libcint cc-pVDZ, NEO/CIAH,
structure-space Davidson, analytic Hessian-diagonal orbital preconditioning,
gradient 1e-3 and energy 1e-7 Eh. All allocate and actually use eight OpenMP
cores within one LLC; BLAS has one thread. Each timing is one native sample,
with saved old data, not a rerun of an old or GitHub executable.

| Scope | Saved run/node | New run/node | Accepted steps | HVP/structure responses | Saved SCF (s) | New SCF (s) | Saved process (s) | New process (s) | Saved/new peak RSS (GiB) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| One step | 1263392/anode50 | 1263446/anode50 | 1 | 6/0 | 7.443532 | 7.055196 | 12.21 | 11.83 | 2.490/2.500 |
| Complete convergence | 1263393/anode21 | 1263447/anode21 | 4 | 38/0 | 31.886365 | 29.082360 | 36.51 | 33.67 | 2.487/2.501 |

The complete process is 7.8% faster; SCF is 8.8% faster. This is a measured
modest endpoint gain, not an order-of-magnitude speedup. HVP counts and every
macro/microstep count are retained (microsteps 10, 33, 19, 3; HVPs 6, 18, 11, 3).
The initial and final energies agree with saved outputs to a few 1e-12 Eh.
New final energy is -399.050612521571 Eh, versus -399.050612521570 Eh saved.
The new run exits zero at the unchanged projected-gradient criterion.
The last energy change is about 5.2e-10 Eh.

RSS is GNU-time process maximum, not Slurm sampling. The increase is about
14 MiB (0.6%); there is no new unbounded array or per-thread carrier. The
explicit CAS12 source-direction workspace is 34,661,760 bytes, within the
unchanged 64 MiB source-workspace bound. Required output tensors and the
molecule-static AO graph are outside that internal workspace accounting.

## Attribution in matched first-step phase captures

Both phase jobs use anode50 and the same eight-core affinity convention.
1263412 captures the current code at entry; 1263453 captures the revised code.
Neither is a native production speed measurement or an old remote rerun.

| Nested phase | Saved calls/time (s) | New calls/time (s) |
|---|---:|---:|
| Orbital HVP | 6 / 2.231320 | 6 / 1.585852 |
| Selected/integral direction | 6 / 1.294953 | 6 / 0.641003 |
| Source sigma | 22 / 1.370182 | 22 / 0.736147 |
| Exterior left transforms | 1444 / 0.243662 | 436 / 0.063829 |
| Separate generator adjoints | 14 / 0.316732 | 8 / 0.222163 |
| Active backward | 6 / 0.289744 | 6 / 0.295129 |
| AO response | 6 / 0.530328 | 6 / 0.532388 |
| Analytic Hessian diagonal | 1 / 2.254031 | 1 / 2.247896 |

Selected-direction time drops 50.5%; HVP time drops 28.9%. The six removed
separate generator calls and 1008 removed exterior calls directly identify
the intended eliminated work. The unchanged expensive AO and diagonal stages
explain why the endpoint reduction is smaller. Phase scopes overlap and cannot
be summed to derive a process-wide percentage.

## Correctness checks and remaining work

Build/test 1263445 passes 7/7 selected CTest entries after relinking their
executables. New scalar tests cover unequal spins, both common exchange
parities, nonunit CI coefficients, nonunit state weights, fixed/moving overlap,
and complete adjoint directions. Their generator oracle is the independent
ordinary full-carrier sigma formula, not the optimized source method. The
product oracle may leave unused lower-triangular generator entries unspecified;
the upper-triangular Cholesky covector and complete physical pullback agree.

CAS12 one-column independent full-product reference errors are 6.36e-15 for
the accepted gradient, 1.30e-18 for H/M images and 1.08e-17 for the complete
directional adjoint. It retains zero accepted/directional product doubles.
The nonclosed three-column product path has directional-gradient finite-
difference error 6.98e-11. These are algebraic/derivative checks, not extra
production certifications. The theory PDF was rebuilt (27 pages), extracted
and the new section visually checked.

The next priority is the **analytic-diagonal AO density contraction**, rather
than another single-column CI micro-optimization. Its needed outputs are
projected density diagonals, not all AO response matrices. Current code has
both canonical packed-scatter and row-owned projection orders; the measured
52-channel path still constructs full AO channel responses, including active
Coulomb responses already available from accepted integrals. Before changing
that kernel, derive exactly which J/K projections are consumed, price their
memory traffic as well as ERI visits, and retain the same signed analytic
diagonal. No such preconditioner or AO-kernel change is bundled into this round.

## Evidence and source scope

`single12_source_jet.csv` and `single12_source_jet_phases.csv` are generated
from logs. `single12_source_jet_evidence.tar.gz` retains native input/output,
hashes, affinity, phase logs, build/tests/oracles, denied profiler stderr, the
scoped patch and Slurm states. Prototype 1263419 had a missing nullptr test
argument; 1263421 used an overstrong full-generator comparison for unconsumed
lower entries; 1263436 passed all CTests but used an incorrect CAS12 runner
argument. Each was corrected before the two production measurements. Dependent
pending jobs were cancelled rather than left unsatisfiable. No old baseline
was rerun. Unrelated RI and gitignore changes were neither transferred nor
included in this round.
