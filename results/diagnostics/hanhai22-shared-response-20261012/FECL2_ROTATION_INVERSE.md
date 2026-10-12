# FeCl2: inter-orbital analytic-diagonal mixing

## Evidence

Full-AO/complete-space analytic-diagonal preconditioning now retains source
mixing before taking the signed inverse. The AH equation, HVP, coordinate
space, additive update, keyframe rules, microstep budgets, and all numerical
thresholds are unchanged. The independent positive block used by block-LBFGS
is unchanged. The factorized construction and its applicability conditions
are derived in `article/neo_rotation_preconditioner.md`.

The local frame inverse uses the existing raw chart Gram factor and is
residual-certified. The signed local eigensystem is released when replaced;
the shared persistent inverse is not duplicated per thread. No full orbital
Hessian or reduced-coordinate matrix is allocated.

Hanhai22 job 1248739, 32 bound CPUs, exact integrals, original FeCl2 input:

| Variant | Steps | HVP | Final energy / Eh | Wall / s | Peak RSS / MiB |
|---|---:|---:|---:|---:|---:|
| Existing overlap-only path, 1248721 | 24 | 544 | -2181.624464442866 | 17.49 | 538.0 |
| Coupled frame preconditioner, 1248739 | 18 | 353 | -2181.624464790287 | 12.02 | 517.9 |

Both converged using the same criteria. HVP count decreased by 35.1%; native
time decreased by 31.3%. Jobs ran on different nodes (anode45/anode18), both
AMD EPYC 9654, so the timing is a single-run cross-node comparison, not a
controlled repeated same-node benchmark. No old version was rerun.

The new endpoint differs from the archived PySCF reference
-2181.6244647028957 Eh by -8.7391e-8 Eh. This is not evidence of a global
minimum or of identical trajectories across implementations.

All 353 Ritz queries had Euclidean retained-space residual divided by full
Euclidean residual below 5.999e-11. Thus solving the projected AH eigenproblem
more accurately does not address these observed complementary residuals.
For macros 15–18, the orbital component accounts for 92.94–98.38% of the
squared **Euclidean coordinate** residual. This is not a decomposition of
the physical dual norm. Tail exits remain `micro_budget` without component
clipping; this improvement does not finish the convergence-efficiency goal.

`test_ciah`, `test_coupled_structure`, `test_orthogonal_direct_ci`, and
`test_sparse_orbital_quotient` passed on the compute node. The last test now
constructs the full frame pullback independently, including inter-orbital
mixing, and compares its inverse to the production factorization. Constrained
HAO tests retain their existing reference.

## Porphine comparison

Final build and unit checks: 1248744. Frozen executable copied into each run.
NEO original-guess runs: 1248745_0 (singlet), 1248745_2 (triplet).
PySCF default mc1step/AH original-guess runs: 1248772_1 (singlet),
1248772_3 (triplet). These tasks are independent and concurrent.

The first PySCF comparison tasks 1248745_1/3 were explicitly cancelled to
test an Angstrom-conversion hypothesis. This hypothesis was **disproved**:
replacement tasks report zero coordinate change relative to the checkpoint,
and normalized AO overlap agrees with the CPP export to 1.055e-15. The
checkpoint already used the CPP loader's 0.529177249 Angstrom/bohr constant.
Different library defaults were therefore not causing this discrepancy.
Replacement tasks keep the explicit same-input Bohr geometry, basis,
Cartesian convention, spin, active space and original orbital coefficients;
no ECP is added.

Initial energies remain -2244.7047401600 (CPP singlet) versus
-2244.70474185271 (PySCF), and -2244.6118593777 (CPP triplet) versus
-2244.61186107004 (PySCF). The common difference is about 1.69e-6 Eh.
Its origin is not established. In particular, scalar ERI screening,
inactive-reference contraction and CI accuracy must not be conflated with
optimization error. Do not rerun fixed-point molecular energy comparisons
to present this hypothesis as resolved.

Both are CAS(8,11), 37 atoms, 469 Cartesian cc-pVDZ AOs, 89 inactive orbitals,
100 stored occupied/active columns. No optimized restart is used. Final
porphine energies, convergence efficiency and RSS are pending; the goal
remains active.
