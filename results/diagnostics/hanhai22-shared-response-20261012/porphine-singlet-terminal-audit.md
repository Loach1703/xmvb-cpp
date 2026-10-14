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

## Corrected original-guess endpoint and orbital-path audit

Corrected singlet job **1253652** completed with 86 accepted macrosteps,
1333 orbital HVPs, energy **-2244.740725398494 Eh**, physical gradient norm
`9.52249589e-4`, and final energy change `-4.6988e-8 Eh`. Wall time was
1:10:46 and peak process RSS 39.5197 GiB. It satisfies the unchanged stopping
contract, but remains **0.344885506 mEh above** the original-guess PySCF
reference. The screening repair did not remove the different-endpoint issue.

Historical files did not retain per-macro orbitals. Array **1253884** therefore
records only the first 1/4 corrected NEO macros and the first 4 PySCF default
AH macros, with the original inputs and thresholds unchanged. The C++ runs
use the same frozen executable as 1253652, the existing scalar event trace,
and a final orbital-table export; no full H/S snapshots or additional
finite-difference audits are enabled. Native exit 2 is the deliberate cap.
PySCF snapshots are passive callbacks; approximate microstep densities are
explicitly distinguished from fully solved accepted macro densities.

Geometric analysis jobs **1253900** and **1253912** evaluate the AO overlap
only and read saved orbitals/densities. They do not solve another CI problem,
evaluate fixed-point energies, or change either optimizer. The active-space
comparisons are invariant to orbital signs, ordering and internal basis changes.

| Diagnostic | Macro 1 | Macro 4 |
|---|---:|---:|
| Largest active-space principal angle, NEO vs default AH | 21.9444 degrees | 89.9569 degrees |
| Largest core-space principal angle | 2.0967 degrees | 2.7250 degrees |
| Physical 1RDM difference, orthonormal AO Frobenius norm | 0.12442 | 0.16594 |

The NEO 1RDM for this comparison is reconstructed from existing natural
orbital output, rounded to six decimals. Its NO Gram error is at most
`3.34e-6` and electron trace error at most `1.39e-5`; these are not
derivative-accuracy measurements. The observed density differences are
well above that output precision.

At macro 4, the nearly orthogonal pair of active principal directions has
occupations approximately `0.000159` (NEO) and `0.000582` (PySCF). Projection
onto the complete Cartesian Fe d AO subspace is respectively `0.986739` and
`6.03e-7`. Thus a weakly occupied Fe-d-like active direction is replaced by
a non-Fe-d direction on the default AH path. These geometric projector
weights are not atomic charges or pure atomic 3d/4d occupations. The data do
not show loss of active-space dimension or an internal orbital gauge error.
They do not alone prove when the two eventual basins became unavoidable.

The first NEO keyframe (after six accumulated microsteps) isolates a model
failure before the first accepted macro:

| Coordinate gradient norm | Predicted | Actual |
|---|---:|---:|
| Orbital | 0.1032800 | 0.1018237 |
| Structure | 0.00859195 | 0.0944827 |

Structure mismatch norm is about `0.0956434`, versus orbital mismatch
`0.0156719`. The accepted macro predicts a decrease of 61.5035 mEh but
actually lowers relaxed energy by only 7.06971 mEh (`rho=0.114948`). It exits
on the existing micro budget, while the predicted gradient is still
0.110992 (source norm 0.164488). These facts locate an early failure of the
finite coupled model's structure prediction. They are not proof of an
incorrect infinitesimal HVP.

Source review confirms an important baseline distinction: the PySCF
reference uses **default mc1step**, not `mc.newton()`. Its AH search is
orbital-only, with approximate CI-density updates between rotations. Our
single-state NEO searches a joint orbital/structure AH space. PySCF's own
coupled `newton_casscf` also retains an origin Hessian across keyframes, so
freezing the Hessian is not by itself a bug. The finite additive CI
normalization pullback in the C++ code includes the normalization factor;
no missing factor has been identified by this source review.

A **one-macro**, original-guess PySCF coupled-Newton control, **1253926_3**,
was submitted to distinguish joint-model path selection from C++-specific
coordinate/update behavior. It is not a new convergence baseline, and no
threshold is varied. Its callback contracts the accepted CI vector for its
1RDM because `newton_casscf`'s macro-local `casdm1` remains stale until
post-loop finalization. Do not use that stale local as accepted density.

Scripts: `porphine-path.sbatch`, `porphine-path-pyscf.py`,
`porphine-path-newton.py`, `porphine-path-analysis.py` and its Slurm wrapper.
Numeric evidence: `porphine-path-1253884/path-comparison.json`.

### Matched coupled-Newton control completed

Job **1253926_3** completed its one diagnostic macro in 130 seconds including
setup/output. Analysis job **1253932** reads its saved orbital frame and CI
1RDM. No energy/gradient reevaluation is performed by the geometric analysis.

| First-macro method | Accepted energy / Eh | Largest active-space angle from initial | HVP/actions recorded |
|---|---:|---:|---:|
| Corrected C++ joint CIAH/NEO | -2244.711811541401 | 19.6137 degrees | 11 orbital HVP |
| PySCF default orbital AH | -2244.722012550095 | 41.1563 degrees | JK count is not an HVP count |
| PySCF coupled Newton | -2244.720783342431 | 37.6230 degrees | 16 coupled HVP |

PySCF coupled Newton's first active subspace is 6.9947 degrees from default
AH, but 18.0391 degrees from C++ NEO. Thus simply attributing the C++
branch to the use of a joint orbital/structure model is insufficient.
PySCF coupled Newton and the C++ code share their headline defaults
(0.03 component cap, 10 base microsteps, AH start cycle 3, keyframe interval 5).
No such parameter was adjusted in this audit.

Remaining source-level differences are the chosen generalized AH metric,
the finite raw-coefficient orbital update, the pulled-back carrier CI
preconditioner and coarse enrichment versus PySCF's component diagonal, and
exact versus approximate keyframe gradients. They can change early finite-subspace
steps. In particular, the C++ generalized metric is a block-diagonal orbital
subspace metric plus structure overlap, not a demonstrated full wavefunction
norm. Nonzero CI representation connection in the Hessian does not by itself
make this positive step-length metric invalid; no missing-cross-metric bug
has been proved.

The diagnosis establishes early active-space branching and concentrated
structure-model prediction error. It does **not** yet causally isolate one
of these remaining differences as the reason for the higher final endpoint.
The next discriminating check is covariance of the first joint AH direction,
metric and preconditioner between the two coordinate representations, not
another endpoint energy replay or a threshold sweep.

### First-keyframe finite-response audit (2026-10-14)

Jobs **1253979** and **1254011** run only the first original-guess singlet
macro, with opt-in diagnostic evaluations. No optimizer parameter, HVP,
convergence criterion or finite update is changed. Both `test_ciah` and
`test_coupled_structure` pass in the Slurm allocation. Timings include the
extra diagnostic evaluations and are not performance results.

The live branch reports `carrier_diagonal=1`. The accepted CI coordinate
gradient norm is `3.71962e-8` (physical dual norm `1.94037e-8`). Thus neither
the raw-VB-diagonal hypothesis nor forcing that tiny origin CI residual to
zero explains the keyframe's approximately `0.09564` structure mismatch.

Let the actual accumulated displacement be $(p,q)$, the origin structure
gradient be $g_{s,0}$, and the two structure components of the origin Hessian
image be $L_o=H_{so}p$ and $L_c=H_{ss}q$. All gradients below are pulled back
to the same immutable origin coordinates. Define

$$
e_o(t)=g_s(tp,0)-g_{s,0}-tL_o,\qquad
e_c(t)=g_s(0,tq)-g_{s,0}-tL_c,
$$

$$
e_m(t)=g_s(tp,tq)-g_s(tp,0)-g_s(0,tq)+g_{s,0}.
$$

The true total remainder is $e_o+e_c+e_m$. That reconstruction is an
algebraic identity, not an independent derivative certification.

| Fraction of the actual displacement | Pure orbital remainder | Pure CI remainder | Mixed remainder | Total remainder |
|---|---:|---:|---:|---:|
| 1/8 | 0.00228273 | 3.97558e-7 | 0.000951480 | 0.00244464 |
| 1/4 | 0.00886515 | 3.18022e-6 | 0.00384984 | 0.00952093 |
| 1/2 | 0.0317566 | 2.54340e-5 | 0.0153382 | 0.0342812 |
| 1 | 0.0889979 | 0.000203224 | 0.0555502 | 0.0956424 |

These are Euclidean structure-covector norms in the C++ chart, not
coordinate-invariant norms shared with PySCF. Physical dual norms show the
same dominance: at the full displacement, the three remainders are
`0.0673105`, `0.000151375` and `0.0416982`.

Two central-difference probes on the actual joint displacement give
structure-HVP errors `1.81120e-8` and `1.68758e-9` at offsets `1e-3` and
`3e-4`. An independent pure-orbital probe at `3e-4` verifies $H_{so}p$ to
`2.99195e-9`. Their image norms are `0.00859195` and `0.0996305`, respectively.
This verifies the linearization along the tested paths, not every Hessian
direction or the finite-step model.

There is strong cancellation in the linear coupled prediction:
$\|L_o\|=0.0996305$, $\|L_c\|=0.0953108$, but
$\|L_o+L_c\|=0.00859195$. Their angle is **175.632 degrees**, and the sum of
individual norms is **22.689 times** the resulting norm. The finite orbital
and mixed terms invalidate that cancellation; the small predicted structure
gradient is consequently not evidence of a small true structure gradient.

Source comparison distinguishes a valid AH metric from a covariant one.
The current block-diagonal metric contains $G_s=Z^TMZ$, with
$\delta C=Zq/\sqrt{2}$ for this single state. In PySCF's instantaneous
orthogonal carrier, $y=F(b)C$, and its normalized tangent is

$$
\delta y=\frac{FZq}{\sqrt{2}}
 +(I-yy^T)\,DF(b)[p]C_0.
$$

Thus the exact pullback of the carrier Euclidean norm has a $G_s/2$ CI
block plus representation-connection cross terms. The current metric is
still an admissible positive step metric, but not that strict PySCF
pullback. This difference can change finite AH directions; it has not yet
been established as the dominant cause.

The C++ raw additive path rebuilds core projection, active overlap and
its Cholesky factor before applying alpha/beta exterior transforms to fixed
raw structure coefficients. PySCF coupled Newton updates its orthogonal
MO frame and normalizes additive CI coefficients in that frame. The finite
paths differ. The existence of the exterior transform alone does not prove
it dominates the measured remainder: orbital variation changes both the
carrier map and the orthogonal integrals. The follow-up geometry audit
measures those representation changes without another CI solve.

Job **1254022** measures the geometry along the same pure-orbital path.
Origin active overlap differs from identity by `3.98932e-14`. Its derivative
along the step has Frobenius norm `0.0459643`. At the full displacement,
the overlap change is `0.367171`, of which the nonlinear remainder has norm
`0.365174`; the Cholesky upper factor changes by `0.198854`. The normalized
carrier coefficient change is `0.0191648`, with nonlinear remainder
`0.0190995`. Its unnormalized norm becomes `0.994614`, so a missing scalar
normalization alone cannot explain the large response mismatch. These
geometric norms do not bound Hamiltonian-weighted gradient changes.

Job **1254034** additionally audits the first PySCF coupled-Newton keyframe.
It leaves the actual approximate-gradient return untouched. It applies the
origin HVP once to the saved cumulative displacement and constructs one
exact gradient at the rotated orbitals with the **same fixed CI vector**;
it performs no extra CI solve. The exact gradient here uses PySCF's
instantaneous rotated coordinate frame, not the C++ raw chart.

| PySCF first-keyframe quantity | Orbital norm | CI norm | Joint norm |
|---|---:|---:|---:|
| Linear prediction | 0.0833668 | 0.0279660 | 0.0879325 |
| Returned approximate keyframe gradient | 0.0822709 | 0.0514441 | 0.0970309 |
| Exact fixed-CI keyframe gradient | 0.0825502 | 0.0530228 | 0.0981119 |
| Prediction minus exact gradient | 0.00856218 | 0.0294426 | 0.0306623 |
| Approximation minus exact gradient | 0.00763778 | 0.00176764 | 0.00783966 |

PySCF therefore also has a finite-step structural-model remainder. Its
small approximate/exact CI-gradient difference (`0.00176764`) does not
support an explanation based solely on comparing an exact C++ keyframe
against an inaccurate PySCF CI keyframe. Do not interpret ratios between
C++ raw-covector norms and PySCF orthogonal-coordinate norms as invariant
error ratios. The saved `first-keyframe.json` and `first-keyframe.npz` retain
the actual PySCF quantities for further mapped comparisons.

### Finite carrier counterfactual: a specific leading source isolated

Job **1254043** adds a matrix-free counterfactual to the same first-macro
audit. Let $F_t$ be the current exact raw-structure to orthogonal
unique-spin-product map, $N_t=\|F_t C_0\|$, and $y_t=F_tC_0/N_t$.
Hold the **orthogonal carrier Hamiltonian $h_0$ fixed**, while allowing only
$F_t$ and its normalization to vary. The resulting structure covector is

$$
g_s^{F}(t)=\frac{\sqrt{2}}{N_t}
 Z^TF_t^T(h_0-E_t^{F})y_t,\qquad E_t^{F}=y_t^Th_0y_t.
$$

This is a coordinate-defined counterfactual, not a separate physical
calculation or a replacement VB wavefunction. On this complete singlet
structure sector, it is evaluated through the retained origin action:

$$
\widetilde C_t=M_0^{-1}F_0^Ty_t,\qquad
h_0y_t=F_0M_0^{-1}H_0\widetilde C_t.
$$

Completeness matters: $F_0M_0^{-1}F_0^T$ is identity on the selected spin
sector, not on arbitrary product vectors or arbitrary selected-structure
subspaces. The measured closure residual
$\|F_0\widetilde C_t-y_t\|$ is at most `9.44e-17` on these probes. Shared
occupation-block inverse factors are reused; no dense H/M matrix, extra
AO integral transformation, CI solve or per-thread carrier copy is made.

The carrier-only remainder is
$e_F(t)=g_s^F(t)-g_s^F(0)-t(g_s^F)'(0)$, with the derivative measured at
`+/-3e-4`. The remaining vector is $e_o(t)-e_F(t)$; it includes Hamiltonian
changes **and their cross terms with the carrier map**, so it is not a
pure-integral contribution.

| Step fraction | Carrier-only remainder | Remaining remainder | Actual pure-orbital remainder |
|---|---:|---:|---:|
| 1/8 | 0.00232136 | 0.000370491 | 0.00228273 |
| 1/4 | 0.00903990 | 0.00168414 | 0.00886515 |
| 1/2 | 0.0324772 | 0.00809489 | 0.0317566 |
| 1 | 0.0898271 | 0.0377999 | 0.0889979 |

These norms do not add: the full-step cross inner product is
`-0.000788553`. No percentage attribution is claimed. Nevertheless, holding
$h_0$ fixed still reproduces the large pure-orbital remainder's scale across
all four fractions. This causally isolates finite carrier-map motion as a
leading mechanism in this chart, rather than blaming electronic-integral
changes or an incorrect directional Hessian alone.

The actual implementation independently changes raw orbital coefficients
and raw structure coefficients. Orbital normalization and core projection
change the active Gram matrix; its Cholesky/exterior map then changes the
normalized carrier CI even with raw $C$ held fixed. The infinitesimal
connection is included in the verified HVP, but its large **finite nonlinear
remainder** is not represented by the origin linear coupled-gradient model.
That remainder spoils the strong orbital/CI cancellation described above.
PySCF's orthogonal-frame finite update does not have this changing-overlap
normalization map. This is a finite-chart/CI-transport issue, not evidence
that nonorthogonal orbitals or VB derivatives are intrinsically incorrect.

The next derivation should align finite orbital retraction and CI transport
while preserving the intended VB variational space and their gradient/HVP
pullbacks. It must distinguish a complete spin sector (exact carrier
transport is representable) from an arbitrary selected VB structure space
(such transport need not remain in the selected space). Merely changing the
CI preconditioner or AH metric does not remove the measured finite map
nonlinearity; no parameter adjustment or production solver change was made.

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
