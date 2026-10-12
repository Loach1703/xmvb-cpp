# FeCl2 frame-step consistency and rejected RHS refresh

All runs use original input orbitals, 32 bound cores on Hanhai22, exact ERIs,
analytic diagonal, and unchanged outer/inner tolerances and step cap 0.03.
Existing baselines were not rerun.

| Variant | Job | Steps | HVP | Final energy / Eh | Wall / s | RSS / MiB |
|---|---|---:|---:|---:|---:|---:|
| Full frame preconditioner | 1248739 | 18 | 353 | -2181.624464790287 | 12.02 | 517.9 |
| Reproject RHS immediately at accepted keyframes | 1248799 | 21 | 390 | -2181.624464855218 | 14.58 | 519.0 |
| Sum source contributions before step-component cap | 1248805 | 16 | 334 | -2181.624464776006 | 11.87 | 516.8 |

The RHS refresh experiment saved some expansion actions at individual
keyframes but increased total work to convergence. Its method, invocation
and experiment-specific unit test were removed; no runtime switch or unused
alternate path remains. Saved native artifacts preserve the evidence.

The frame-step correction retains the actual linear frame-map convention:
each source orbital contributes an outer product to a common external-by-source
rotation matrix; those products must be added before taking the largest entry.
The old max-of-individual-products expression need not equal that entry.
The cap is still 0.03 and the orbital update remains additive in the same chart.
Constrained local charts retain their existing convention. The full-frame
unit reference uses the independently assembled linear map on the total step.

The new endpoint differs by -7.3110e-8 Eh from the archived PySCF reference
-2181.6244647028957 Eh. `test_ciah`, `test_coupled_structure`,
`test_orthogonal_direct_ci` and `test_sparse_orbital_quotient` passed in
1248805. Native times are single runs on different nodes with the same EPYC
9654 CPU model. The small 0.15 s timing change is not a robust speedup claim.
Tail micro-budget exits remain, so efficient convergence is still unresolved.

Next structural check: the diagonal's virtual reference frame currently comes
from the S-metric complement decomposition (`preparer.cpp`), not an effective
one-electron spectral frame. A diagonal curvature model is not invariant
under a general rotation of that complement. Investigate the missing external
row couplings before changing numerical stopping policies.
