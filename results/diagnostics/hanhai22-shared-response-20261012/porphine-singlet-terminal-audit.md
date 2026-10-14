# Iron porphine singlet: terminal-point and integral audit

2026-10-14. Runtime optimizer: master `8e9d8e1e`; Hanhai22, 32 cores,
exact four-center integrals, Cartesian cc-pVDZ, CAS(8,11), 89 inactive
orbitals, 469 AO functions, no ECP. No optimizer or stopping thresholds were
changed in this audit. All quantum-chemical calculations used Slurm compute nodes.

## Energy endpoints

| Calculation | Energy / Eh | Macro steps |
|---|---:|---:|
| NEO before shell-screen fix, original guess, job 1252776 | -2244.740723326940 | 81 |
| PySCF default AH, original guess, job 1248772 | -2244.741070284000 | 128 |
| PySCF default AH, NEO terminal wavefunction, job 1253506 | -2244.740725080374 | 1 additional |

PySCF's restart starts at -2244.7407250192023 Eh and decreases by
6.11717e-8 Eh before declaring convergence under the unchanged thresholds.
The 0.346957 mEh original gap is therefore not removed by handing the NEO
terminal wavefunction to PySCF. The remaining difference between the two
PySCF endpoints is 0.345204 mEh. This is evidence of different approximate
stationary endpoints, not proof that either is a strict local minimum.

The checkpoint active 1RDM is stored in the same active-orbital gauge as its
MO coefficients: default PySCF canonicalization does not rotate active
orbitals when `natorb=False`. Comparing both physical spin-summed 1RDMs in
the same orthonormal AO frame gives a Frobenius distance of 0.6450916.
Thus the different active subspaces are not merely a harmless orbital gauge
change. Their most weakly occupied active natural orbitals have occupations
0.00127957 (NEO endpoint) and 0.00079356 (PySCF reference).

## Coordinate and derivative checks

- The full singlet structure space contains 32670 CSFs. The orbital chart has
  37879 directions: 979 inactive-active, 32841 inactive-virtual and 4059
  active-virtual. Its dimension matches the complete CAS variation space.
  Dimension alone does not establish numerical correctness of its map.
- The audited endpoint replay differs from the reported NEO terminal energy
  by only 2.6e-11 Eh; raw coefficients, raw gradients and exported carrier
  came from the same accepted evaluation.
- A sampled coupled HVP agrees with central differences of coupled keyframe
  gradients; decreasing the displacement from 1e-3 to 3e-4 decreases orbital
  error from 1.95608e-5 to 1.76083e-6 and structure error from 3.03775e-5 to
  2.73397e-6. This validates the sampled direction, not the entire operator.
- After analytic conversion into all 37879 orthogonal rotation components,
  C++ and PySCF gradients differ by norm 3.21848e-4, almost entirely in the
  inactive-virtual sector (3.21802e-4). The orbital representation error in
  this conversion is 2.8e-15.
- Targeted energy differences in this sector approach the C++ derivative,
  not the PySCF derivative. The smallest displacements become contaminated
  by CI/energy resolution and cannot establish high relative derivative
  accuracy. No CI or optimizer tolerance was tightened to hide this issue.

## Where the gradient difference originates

Job 1253529 uses the same AO density for both Fock constructions.

| Quantity | Maximum absolute difference |
|---|---:|
| AO overlap, unit-normalized convention | 1.05e-15 |
| Inactive density matrix | 1.78e-15 |
| AO core Hamiltonian | 1.15e-11 |
| Inactive Fock matrix | 9.28124e-7 |

Projecting four times the inactive Fock difference into inactive-virtual
directions gives norm 3.20525e-4, accounting for nearly all of the full
gradient discrepancy. Thus an AO two-electron operator difference is already
present before the orbital pullback. This does not identify its root cause.

## PySCF screening test

Job 1253552 used installed PySCF 2.14.0. No RI/density fitting or ECP was
enabled. Configuration defaults are `direct_scf_tol=1e-13` and AO2MO Schwarz
cutoff `1e-14`; these are not a guarantee that those branches execute.
The diagnostic actually materializes a full in-core ERI array, so its J/K
contraction does not use the direct-SCF prescreen. AO2MO also has a separate
in-core branch that transforms the available ERI array.

For the AO pair (185,454), zero-based, two fixed-index ERI slabs were evaluated
and contracted directly. This avoids the J/K optimizer's screening and never
forms a full four-index diagnostic tensor. The explicit libcint exponential
screen was set to 1e-16 for this slab calculation; this is not a claim of
mathematically untruncated integrals.

| Fock entry in the PySCF AO convention | Value |
|---|---:|
| PySCF standard construction | 0.0005318013826473075 |
| Direct ERI-slab contraction | 0.0005318013826466130 |
| C++ construction | 0.0005308732587864763 |

The PySCF/direct difference is 6.95e-16, versus a C++/direct difference of
-9.28124e-7. Applying the C++ per-integral 1e-10 value discard rule to these
same slabs predicts a bias of only +3.14583e-10, with the opposite sign.
Therefore neither PySCF's J/K screening nor the C++ individual-integral
discard alone explains this particular Fock discrepancy. Shell-pair
prescreening, integral generation and graph contraction were then isolated in
the subsequent test below; no threshold change was needed.

## Confirmed shell-screening error and repair

The independent scalar probe in job 1253595 reconstructs the full ordered
density contraction from shell quartets, without building a global AO graph.
Using the old raw-shell screening reproduces the original Fock component to
1e-17. Removing only that screening produces the PySCF/direct result to 1e-17.
The graph contraction and four-center evaluator are therefore not the cause
of this component's discrepancy.

The old shell-pair bound was computed before per-AO normalization, whereas
the retained quartets were subsequently normalized. The largest AO factor
in this basis is 7.12064. A bound on raw integrals does not bound the normalized
integrals at the same numerical threshold. In this scalar contraction, the
old screening incorrectly discarded 18624 Coulomb and 5204 exchange terms
whose normalized integral magnitude exceeded 1e-10; their largest magnitudes
were 9.99e-8 and 7.80e-8, respectively.

The repaired bound is the maximum self-Coulomb integral of the *normalized*
AO products in each shell pair. The square root of the product of two such
maxima is a valid Schwarz upper bound on every normalized quartet component.
The implementation reads only the self-diagonal and applies its two AO
normalization factors squared. It does not alter the integral cutoff, the
optimization tolerances, the orbital chart, or the HVP definition. The raw
bound API and misleading historical option name were removed.

Job 1253607 validates the repair:

- The screened scalar Fock entry is 0.00053180169722949141, exactly matching
  the independently value-screened result. No integral above the cutoff is
  missed by shell screening. The remaining +3.14583e-10 discrepancy is the
  already-quantified per-integral discard contribution.
- Unit regressions cover all s/p/d quartets, the Schwarz inequality, invariance
  under contraction-coefficient rescaling, and equality of screened/unscreened
  materialized streams. The new regression passes with 1 and 4 OpenMP threads;
  the existing libcint RI regression also passes.
- The audited terminal energy is -2244.7407249986918 Eh, differing from the
  PySCF energy for that same orbital point by 2.05e-8 Eh instead of 1.69e-6 Eh.
- The sampled coupled HVP still agrees with keyframe gradient differences;
  reducing displacement from 1e-3 to 3e-4 reduces orbital error from
  1.77417e-5 to 1.59705e-6 and structure error from 2.92730e-5 to 2.63456e-6.

Job 1253624 compares all orbital derivatives with the same corrected C++ CI
wavefunction in PySCF. The gradient difference norm drops from 3.21848e-4 to
8.52863e-6 (about 38-fold), with maximum component error 7.58587e-7. It is not
machine-precision agreement: the unchanged per-integral discard remains,
and its full-vector contribution has not been separately reconstructed.
The corrected state is a pure singlet (spin square about 3e-15).

This fixes an integral-generation error, not the different-basin issue. It
does not establish that the corrected optimizer reaches the lower PySCF
endpoint from the original guess. Job 1253625 runs small molecular regressions
and a four-macro original-guess porphine diagnostic to investigate that next.

| Corrected native run | Macros | Orbital HVP | Final energy / Eh | Wall / s | Peak RSS / GiB |
|---|---:|---:|---:|---:|---:|
| F2 OEO | 3 | 12 | -198.689800055476 | 0.17 | 0.018 |
| FeCl2 | 13 | 270 | -2181.624464864732 | 5.49 | 0.497 |
| MnF2 | 11 | 158 | -1348.893353292947 | 4.15 | 0.531 |
| Porphine singlet, four-step diagnostic only | 4 | 41 | -2244.717160198732 | 245.98 | 39.519 |

The three small systems converged. The porphine diagnostic deliberately stops
at the iteration limit (native exit 2), not a numerical failure. Its first
three steps remain close to the old trajectory, while the fourth differs;
both use 41 HVPs over those four steps. This is not evidence of faster or
lower-energy convergence. Job 1253652 now runs the corrected singlet to
convergence from the original guess, not from a saved endpoint. Only the
integral screening correction differs from the previous original-guess run.

## Molcas third-party calculation

Job 1253591 runs native RASSCF/QUNE for the singlet from the same original GUS
frame, CAS(8,11), 89 inactive orbitals, 469 Cartesian AO functions, no ECP,
no RI/Cholesky. QUNE is Molcas's native super-CI quasi-Newton update, not the
same solver as C++ NEO or PySCF AH.

The old Molcas converter failed because the native Fe cc-pVDZ generalized
contractions mix in an uncontracted diffuse primitive; a permutation/scale-only
map cannot represent this change. The new input writes the actual reference
contraction columns into a custom nonrelativistic basis library, rather than
silently using a different AO convention or an endpoint guess. ONEONLY setup
checks conversion before the costly four-center integral generation.
Normalized AO overlap error is 1.37e-13 and orbital Gram error 1.63e-12.
The job has entered conventional integral generation; no final Molcas energy
is available yet.

## Artifacts

- `porphine-terminal-1253506/cpp/audit.out`: endpoint chart and coupled HVP audit.
- `porphine-terminal-1253506/pyscf/gradient.json`: cross-program direction check.
- `porphine-terminal-1253506/pyscf/output.out`: authoritative macro trajectory.
  The callback records both micro and macro events; its first `energy_change`
  in `result.json` is stale and must not be interpreted as a macro decrease.
- `porphine-gradient-map-1253511/comparison.json`: all-component gradient mapping.
- `porphine-ao-adjoint-1253529/ao-comparison.json`: AO-level isolation.
- `porphine-pyscf-screen-1253552/screening.json`: runtime flags, ERI slab control,
  natural occupations and physical-density comparison.

Jobs 1253518 and 1253549 failed in Python setup/reporting, not the native
optimizer. Job 1253527 completed the C++ finite differences and exports, then
failed its Python reporting phase after a module-path switch; its C++ output
remains usable. The later AO/slab checks used those completed exports.
