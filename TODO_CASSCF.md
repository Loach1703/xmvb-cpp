# Orthogonal CASSCF

## Scope and ownership

Reuse one CIAH solver with separate physical-model adapters. Do not clone the
VB macro driver, construct VB structures for CASSCF, or introduce compatibility
headers. The initial scope is real, state-specific, full-CAS, exact integrals.

- `mcscf/ci`: physical active integrals, spin strings, sigma/adjoint and spin.
- `linalg`: eigensolvers, eigen response, metric solves and Krylov algebra.
- `parallel`: shared-buffer OpenMP execution and work partitioning.
- `integrals`: common AO operators, exact/RI active transforms, integral
  directions and adjoint/HVP kernels, without VB input dependencies.
- `mcscf/optimization`: shared CIAH and norm-extended subspace solvers and algebraic
  preconditioners; vector actions and dual norms are supplied by each model.
- `mcscf/states`: one implementation of state-weight normalization.
- `casscf/integrals`: inactive energy and active-Hamiltonian construction.
- `casscf/driver`: orthogonal-model state and CI solve.
- `vbscf`: nonorthogonal charts, selected structures and their transport.

The common and CASSCF targets must not include VBSCF headers or link to the VB
library. CMake checks this dependency direction. An active integral provider
must not require a full four-index MO integral array.

## Milestones

1. Shared direct-CI and independent pure-spin CASCI.
   - Implemented: common sigma, physical integral data, matrix-free spin
     projection, inactive/nuclear constant and existing core Davidson.
   - Complete: Hanhai22 job `1257606`, 15/15 regressions passed, including
     original VB direct-CI, CIAH, transport and F2 optimizer references.
   - Independent PySCF-integral CASCI references: F2 CAS(2,2) singlet,
     O2 CAS(8,6) triplet, H4 CAS(4,4) singlet, all STO-3G. Energy differences
     were below 3e-13 Eh. These are fixed-orbital tests, not CASSCF convergence
     or performance results.
   - Pure-spin projection/idempotence, multiple roots, inactive/nuclear energy,
     alpha/beta imbalance and pre-allocation resource rejection are covered.
2. Orthogonal rotations, analytic gradient, fixed-CI HVP and analytic diagonal.
   - Integral prerequisites extracted: exact/RI AO actions, forward active
     transforms, matrix pullbacks and integral-response actions now live in
     `xmvb_integrals`; VB preparation is an adapter. The AO-to-CASCI entry
     consumes these same kernels, with an independent rotated-MO reference.
   - This does not yet implement orthogonal orbital-coordinate derivatives or
     a complete CASSCF macro driver.
   - Use inactive/active/virtual cross-sector rotations only for full CAS.
   - Verify directional derivatives, HVPs and bilinear symmetry.
3. Joint orbital/CI action and CI tangent.
   - Preserve target spin in guesses and every correction direction.
   - No additional relaxed-Hessian response solve inside the joint action.
4. Shared CIAH macro policy, with VB and orthogonal model callbacks.
   - Numerical solvers and their tests are now independent of VBSCF. The VB
     macrostep driver still owns charts, CI transport and physical actions;
     extracting that policy remains separate work, not a second solver copy.
   - First extract policy from the working VB driver; do not copy it.
   - Include finite-update derivative pullbacks consistently.
5. Same-guess validation: F2, 241, FeCl2, then Fe-porphine.
   - Compute on Hanhai22 Slurm nodes; retain original input/binary hashes.
   - Report energy, macros, actual action counts, wall time and peak RSS.
6. Equal-weight state averages on the same interface.
   - Check state-space gauge and spin, not just averaged energies.

## Coordinate contract for milestones 2–4

Two valid conventions must not be mixed. For a fixed origin chart,

$$
B(x)=B_0\exp(\mathcal Lx),\qquad B_0^TSB_0=I,
$$

where `L` packs the independent parameters into an antisymmetric generator.
If `G` is the current right-local orbital gradient, the fixed-chart gradient is

$$
g_x=\mathcal L^*\mathcal J_\kappa^*[G],\qquad
\mathcal J_\kappa[A]=\int_0^1e^{-s\kappa}Ae^{s\kappa}\,ds.
$$

For a fixed CI tangent `Q`, the normalized coefficient map is

$$
z=C_0+Qq,\qquad C(q)=z/\|z\|,\qquad
DC=(I-CC^T)Q/\|z\|.
$$

The complete retraction `R(x,q)=(B(x),C(q))` requires

$$
\nabla^2(E\circ R)[v]
=J_R^T\nabla^2E[J_Rv]+(DJ_R[v])^T\nabla E.
$$

These are update-map terms, not additional physical Hessian blocks. Current
nonorthogonal Gram-driven CI transport becomes an ordinary CI increment when
the active orbitals remain orthonormal; it is not a new CASSCF transport term.

Alternatively, a moving-chart CIAH commits `U <- U exp(delta)` and resets the
segment displacement at a keyframe, as PySCF does. Its local gradient must not
be mistaken for a derivative with respect to one fixed cumulative logarithm.
Reusing a frozen curvature in either convention is an explicit model
approximation, not an exact updated Hessian. Select and test the orthogonal
adapter's convention before milestone 2; do not mix them during implementation.

## Candidate transferable improvement

Retain block-separated samples `(p,0)` and `(0,q)` and their existing coupled
images. For the shifted operator `K=H-omega G`, `E=Z^TKZ`, and
`Q_c=ZE^{-1}Z^T`, the current balancing inverse is

$$
P=Q_c+(I-Q_cK)P_0(I-KQ_c),\qquad PKZ=Z.
$$

This algebra is general two-level preconditioning. The specific candidate is
the block-separated coarse space obtained without additional HVPs, retaining
orbital–CI coupling. It survives the orthogonal limit. Do not equate it with a
new balancing identity or claim a CASSCF speedup before an isolated comparison.
