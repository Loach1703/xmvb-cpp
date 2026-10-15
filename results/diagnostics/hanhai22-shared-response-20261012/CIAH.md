# NEO with CIAH microiterations

## Scope

The single-state NEO backend now uses the coupled CIAH workflow described by
PySCF's `mcscf/newton_casscf.py` and `soscf/ciah.py`. CIAH is the iterative
solution/update strategy for the augmented-Hessian model, not a replacement
by a first-order optimizer. The nonorthogonal VB orbital quotient, physical
metric, analytic coupled HVP and additive raw-coefficient retraction remain
the production VBSCF operators.

Sources consulted:
- https://pyscf.org/_modules/pyscf/mcscf/newton_casscf.html
- https://pyscf.org/_modules/pyscf/soscf/ciah.html

## Workflow

1. Build one immutable coupled orbital/structure model at the accepted point.
2. Stream Davidson directions and retain their Hessian and metric images.
   Refresh only the augmented gradient coupling after a microstep or keyframe.
3. Select the first ascending augmented root with reference weight above 0.1,
   rather than separately certifying the lowest eigenvalue of the Hessian.
4. Scale the orbital/structure microstep and its cached Hessian image together.
   Update the predicted gradient without another HVP.
5. Evaluate fixed-structure, nonlinear keyframes when the inherited CIAH policy
   requests them. Pull both gradients back to the same frozen coordinates;
   include the finite structure-coefficient normalization derivative.
6. Roll back the most recent microstep if the gradient prediction fails the
   keyframe guard. Rebuild the model only after a fully relaxed macro endpoint.
7. Retain the last *accepted* orbital and structure direction as the next
   Davidson initial guess. Do not clear it merely because the outer gradient
   has become small. This is subspace initialization, not a secant history.

Initial microstep policies use PySCF newton-CASSCF defaults (component bound
0.03, three-direction early start, startup residual 500, base micro limit 10,
keyframe interval 5, gradient trust factor 3, Davidson space 30). These are
inner policies, not changes to the outer gradient/energy accuracy contract.
The VBSCF energy-acceptance check remains in place. Backtracking scales the
accumulated coupled displacement and reuses its Hessian image.

Automatic single-state preconditioning uses the signed analytic Hessian model.
Explicit one-electron preconditioning retains its local signed curvature;
the analytic physical step map is supplied independently. BLAS stays
single-threaded inside the 32-thread OpenMP run. Large arrays are not copied
per thread, and the Davidson workspace stores at most 30 vector/image columns.

State-averaged NEO retains its existing implementation: finite multistate
coefficient orthogonality must be treated explicitly before porting it to this
single-state keyframe implementation. The existing cost-based admission rule
can also choose orbital-only NEO for large structure spaces. FeCl2 admits the
complete coupled model (3599 orbital plus 419 structure coordinates).

## FeCl2 evidence

Exact libcint integrals, Cartesian def2-TZVP, FULL/OEO, CAS(10,8), quintet,
Davidson. Hanhai22, 32 bound cores, gradient 1e-3 and energy change 1e-7 Eh.
No endpoint stability criterion was added; no old executable was rerun.

| Version/job | Macro steps | Orbital HVP | Energy (Eh) | Native wall (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|
| Prior correct NEO, 1247396 (saved data) | 17 | 5804 | -2181.624464422983 | 127.49 | 492.67 |
| Initial CIAH integration, 1247507 | 26 | 606 | -2181.624462008713 | 15.94 | 414.93 |
| Complete direction recycling, 1247527 | 22 | 483 | -2181.624464384744 | 14.75 | 428.77 |
| Final source-matched CIAH, 1247545 | 21 | 457 | -2181.624464192990 | 14.35 | 401.89 |

The first CIAH integration is not the accepted numerical result. It inherited
the old small-gradient direction reset and discarded the structure component
of the next initial guess. Fixing these integration differences restores
agreement with the prior correct NEO endpoint (difference 3.82e-8 Eh).
The PySCF reference is -2181.6244647028957 Eh; the difference is 3.18e-7 Eh,
comparable to the already accepted 2.80e-7 Eh of the prior correct run.
The final build additionally matches PySCF's continued residual expansion
after a yielded converged candidate, its eigenvalue-change history across
microsteps, and the inner gradient tolerance passed by the outer driver.
Its final difference from PySCF is 5.10e-7 Eh; from the prior correct NEO it
is 2.30e-7 Eh. This is sub-microhartree agreement, not equality to every
printed decimal. Outer tolerances are unchanged. HVP count falls by 12.7x;
native wall time falls by 8.9x in these single samples.

Final-build default regressions, array 1247549:

| System | Macro steps | Orbital HVP | Energy (Eh) | Native wall (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|
| F2 | 3 | 24 | -198.751155830525 | 0.28 | 17.39 |
| 241 | 3 | 26 | -230.720590392853 | 2.43 | 821.55 |
| MnF2 | 11 | 169 | -1348.893353223711 | 5.49 | 439.59 |

The extra FeCl2 **explicit one-electron** experiment in 1247549_3 was stopped
after more than 100 macro iterations on the old high-energy branch near
-2181.6176388 Eh. It is not a successful regression and is not used as timing
evidence. The explicit one-electron initialization bug found in 1247534_3
(unconstructed inactive projector) was fixed before this final test. The
result shows that the preconditioner can change finite-subspace coverage;
the CIAH reference-weight rule is not a proof of global saddle escape.
Do not claim uniform robustness of all selectable preconditioners.

Job 1247527 used anode18, while 1247507 used anode04. These are the same
Hanhai22 CPU partition, but timings are single native samples, not precise
hardware-matched speed factors. The dominant reduction is directly measured:
5804 to 483 orbital HVPs. The structure-action counter now counts direct
coupled CI images; it is not equivalent to the old iterative response-solve
counter.

Raw results include `output.out`, `neo.tsv`, GNU `resource.time`, input and
binary hashes, node/CPU metadata, build flags, module/thread settings and unit
test logs. `test_ciah`, `test_neo_response_solver` and `test_coupled_structure`
passed before the chemistry run. Job 1247503 stopped during build target
discovery after CMake regeneration; it did not run chemistry.
