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

### Projected finite CI transport: derivation and first-keyframe validation

The theory is recorded separately in `article/finite_ci_transport.md`.
It preserves the existing nonorthogonal raw orbital update and its origin
first derivative. It retains the origin carrier connection in the target,
then projects that target into the **current selected VB image**:

$$
v=F_0(C_0+Zq/\sqrt2)+DF_0[p]C_0,\qquad
M_p C_*=F_p^T v.
$$

For arbitrary selection, the projection residual is generally nonzero.
The derivative includes $(dF_p)^T(v-F_p C_*)$; it is not dropped. Finite
keyframe gradients require the full coefficient-transport chain, including
its orbital term. At the origin, an approximate CI state also requires a
residual-weighted Hessian correction. These terms are derived explicitly;
the corrected finite update must not be enabled with the old fixed-C
keyframe gradient alone.

Job **1254085** builds the diagnostic and runs three unit executables on
Hanhai22 compute node anode01, then audits the original singlet's first
keyframe. The optimizer continues to consume its original candidate;
the transported values are **not** returned to its AH workspace. The AH
metric, preconditioner, origin HVP and outer thresholds are unchanged.
This is not a convergence or performance benchmark.

| Same first-keyframe orbital displacement | Original raw CI update | Projected finite transport |
|---|---:|---:|
| Actual structure gradient norm, origin q coordinates | 0.0944827 | 0.0322300 |
| Error against the same origin structure HVP prediction | 0.0956424 | 0.0267007 |
| Fixed-CI total energy / Eh | -2244.711209544314 | -2244.712506527170 |

The mismatch is reduced by **3.5820 times**. Both gradients use the same
origin q coordinates; the transported gradient includes its exact q-chain
pullback. The energy row compares prescribed CI states at the identical
orbital keyframe, not independently relaxed CI energies or final endpoints.
The transported coefficient correction norm is `0.0260300`. The complete
sector's measured carrier projection residual is `5.88710e-17`, normal
residual `1.12284e-16`, and relative metric-solve residual `1.04131e-16`.

Derivative checks along the actual accumulated direction give:

| Check | Result |
|---|---:|
| Structure HVP error, offset 1e-3 | 1.65387e-8 |
| Structure HVP error, offset 3e-4 | 1.55013e-9 |
| Included CI-residual correction to structure HVP | 1.02289e-9 |
| Origin coefficient first-derivative error, offset 3e-4 | 4.51005e-10 |
| Finite joint energy derivative: analytic | -0.00729314316946 |
| Finite joint energy derivative: central difference | -0.00729317738054 |
| Joint derivative absolute discrepancy | 3.42111e-8 |
| Transport coefficient directional-derivative discrepancy | 2.85204e-10 |
| Separate finite CI energy derivative discrepancy | 8.20407e-10 |

The finite joint derivative contains both the fixed-C orbital derivative
and $g_C^T DC_*$. The CI derivative uses $L^TF_0^TF_pM_p^{-1}g_C$, not the
old fixed-C projection. Only the structure component of the molecular HVP
is certified by this probe; a production full-vector transported HVP still
needs the orbital residual-adjoint terms described in the theory.

`test_ci_transport` also uses a rectangular synthetic map with a deliberately
nonstationary CI state to verify the Hessian correction, and a genuine
three-orbital exterior map with **three selected doubly occupied structures**
(two electrons in each spin sector). Its selected carrier is 9 by 3, and
is not invariant under the exterior transform. Maximum unit-check error is
`4.68039e-8`; the reported nonclosure residual is `0.0140568`. Normal-equation
stationarity and the finite CI gradient pass despite nonclosure. Thus the
closed molecular example is not being used as a proof for arbitrary subsets.

All three executables (`test_ci_transport`, `test_ciah`,
`test_coupled_structure`) pass. The native molecular application exits 2
at its intentional one-macro cap; Slurm wrapper completes successfully.
It performs nine additional fixed-C gradient evaluations, costing 29.4695 s
in total. GNU process wall time is 136.37 s and peak RSS 41,438,504 KiB
(39.5188 GiB), including those audits. These are diagnostic resources, not
the cost of an integrated production transport policy. No AO ERI array is
copied into the trial action or per-thread workspaces.

The short test supports the finite-transport mechanism and a substantially
better first structural prediction. It does not yet establish a reduced
macro count or convergence to the PySCF endpoint. Before any production
integration, complete the transported orbital-gradient adjoint and the
CI-residual full HVP terms; arbitrary subsets require their selected-space
overlap solve, not the complete-sector factor inverse used by this admitted
molecular prototype.

Artifacts: `porphine-ci-transport-1254085/run.err`, `neo.tsv`,
`resource.time`, source/binary hashes and unit outputs. The bounded wrapper
is `porphine-ci-transport.sbatch`; diagnostic mode is the existing audit
variable with value `transport`, not a new optimizer parameter.

## Production finite transport and full pullback (2026-10-14)

The prototype has now been replaced by `optimization/coupled/CiTransport`.
It uses actual selected-space cross-Gram actions and normal equations, not
an assumption that the selected structure space is closed. Complete ambient
spin-string carriers use exterior-transform actions; incomplete carriers use
directed, tiled unique-string overlap contractions without full-space padding.
The actual trial transforms reuse immutable complete-sector metric factors.
No origin factor is rebound and no AO ERI array is copied per thread.

Finite keyframe gradients include normalization, the inactive projector,
active-overlap Cholesky derivatives, and both sides of the normal equation.
The origin HVP retains the nonstationary CI residual curvature. The old
diagnostic-only transport source has been removed. Coordinates, AH policy,
preconditioner, screening, and convergence thresholds are unchanged.

Build 1254176 passed all eight tests, including independent cross-Gram
adjoints, raw-overlap/Cholesky pullbacks, and an arbitrary three-dimensional
selected-space transport test with a deliberately nonstationary CI vector.
The latter verifies both orbital/structure gradients and the full origin
curvature, not just coefficient transport.

The first molecular keyframe was audited using the actual accumulated
orbital/structure displacement. Absolute HVP errors below are in the full
joint coordinate vector; the gradient check compares a finite-point energy
directional derivative to the full transported gradient dot displacement.

| Molecule | FD scale | Full HVP error | Orbital HVP error | Structure HVP error | Energy derivative error (Eh) |
|---|---:|---:|---:|---:|---:|
| FeCl2 | 1e-3 | 4.7804e-10 | 3.4203e-10 | 3.3397e-10 | 4.1929e-9 |
| FeCl2 | 3e-4 | 1.0024e-9 | 9.7616e-10 | 2.2770e-10 | 2.5290e-8 |
| Porphine singlet | 1e-3 | 1.6885e-8 | 3.2314e-9 | 1.6573e-8 | 1.1698e-8 |
| Porphine singlet | 3e-4 | 6.8410e-9 | 6.6722e-9 | 1.5103e-9 | 1.9679e-8 |

Native 32-thread Hanhai22 small-system results (GNU process wall/RSS;
diagnostic audits disabled):

| Case | Macro steps | Orbital HVP | Final energy (Eh) | Wall (s) | RSS (MiB) |
|---|---:|---:|---:|---:|---:|
| F2-OEO | 3 | 15 | -198.689800050559 | 0.30 | 18.51 |
| FeCl2 | 26 | 455 | -2181.624464867798 | 11.69 | 495.82 |
| MnF2 | 11 | 157 | -1348.893353292904 | 3.87 | 499.16 |
| 241 | 3 | 20 | -230.720590533847 | 2.01 | 810.52 |

Existing screening-corrected native runs 1253625 are reused, not rerun:
F2-OEO 3 macros/12 HVP/0.17 s, FeCl2 13/270/5.49 s, MnF2 11/158/4.15 s.
Final energy differences are respectively +4.917e-9, -3.066e-9, and
+4.3e-11 Eh. Thus energy parity holds, but this transport is **not a universal
convergence improvement**: FeCl2 takes more steps. Independent integration
reviews found no origin gradient scaling or duplicated pullback terms;
the molecular FD checks support the implemented derivatives. Step-count
regression must not be hidden by changing tolerances.

For porphine singlet, the first macro starts at -2244.7047418329 Eh in both
runs. Existing run 1253652 reaches -2244.7118115414 Eh with 11 HVP;
integrated transport reaches -2244.7197611390 Eh with the same 11 HVP.
Projected infinity gradient changes from 1.20833749e-2 to 1.05294761e-2.
This is a 2.12x larger first-macro energy decrease, not a convergence proof.
The audited one-macro run costs 137.31 s and 39.521 GiB RSS, including eight
additional gradient evaluations; it is not a production timing comparison.
Its application exit 2 is the intentional diagnostic iteration cap.

After short validation, native original-guess convergence jobs were submitted
as array 1254189_0 (singlet) and 1254189_1 (triplet), 32 cores each. No audit
evaluations or shortened optimization cap are enabled. Both initially queue
on scheduler priority. Their final energy, macro/HVP count, wall time, and
peak RSS remain pending. The frozen binary SHA256 is
`6c1fa8192000f4e5c2d7aed8879f3710394ef10354d28b6a6268c99e3e0d247c`.

Artifacts: `transport-build-1254176`, `transport-regression-1254179`,
`transport-regression-1254187` (FeCl2 first-macro audit),
`transport-build.sbatch`, `transport-regression.sbatch`, and
`porphine-transport-converge.sbatch`. The theoretical derivation remains in
the local, ignored `article/finite_ci_transport.md`.

## Actual raw orbital angles and radial components (2026-10-14)

This is a read-only geometry diagnostic, not an exponential-update experiment.
The raw update, finite CI transport, AH policy, Hessian, preconditioner, and
accuracy criteria are unchanged. The existing long-run frozen binaries are
untouched. Build 1254268 passes the geometry unit test, CIAH, and selected-space
CI transport. Geometry logging is default-off and enabled only by
`XMVB_NEO_ORBITAL_STEP_AUDIT=1` with scalar trace. No additional electronic
energy/gradient/HVP is evaluated. The AO overlap is borrowed and AO vectors
are reused serially; no per-thread metric matrix is created.

For each occupied raw orbital, the recorded definitions are

$$
n_i=\sqrt{b_i^TSb_i},\quad x_i=b_i/n_i,\quad
d_i=\delta b_i/n_i,\quad \alpha_i=x_i^TSd_i,\quad
v_i=d_i-\alpha_i x_i,\quad \theta_i=\sqrt{v_i^TSv_i}.
$$

The directly normalized raw endpoint has angle

$$
\phi_i=\operatorname{atan2}(\theta_i,1+\alpha_i),\qquad
\frac{\|b_i+\delta b_i\|_S}{n_i}
=\sqrt{(1+\alpha_i)^2+\theta_i^2}.
$$

All angles are radians **before core projection**. They are not active-space
principal angles or the full AH metric norm. `micro_origin` measures the clipped
AH increment at the macro origin; `micro_increment` measures that same increment
at the current accumulated raw coefficients. `keyframe`, `macro_trial`, and
`macro_accepted` measure cumulative displacements from the immutable macro origin.
Trial and keyframe samples may subsequently be rejected or shortened.

Array 1254271 uses 32 cores on Hanhai22. FeCl2 is cheap enough to record to
convergence; porphine states are intentionally capped at two macros. Statistics
below use **accepted macro steps and active orbitals only**. The weighted radial
fraction is $\sqrt{\sum\alpha_i^2/\sum(\alpha_i^2+\theta_i^2)}$, summed over these
samples. It is an AO-normalized displacement diagnostic, not a physical norm.

| Case / scope | Samples | Max abs(alpha) | theta p95 | Max theta | Weighted radial fraction |
|---|---:|---:|---:|---:|---:|
| FeCl2, 17 macros | 136 | 0.150622 | 0.125087 | 0.269591 | 35.294% |
| Porphine singlet, first 2 macros | 22 | 0.565578 | 0.831260 | 1.206744 | 29.172% |
| Porphine triplet, first 2 macros | 22 | 0.977439 | 1.290013 | 4.767531 | 23.172% |

Concrete same-orbital accepted examples (orbital indices are one-based):

| Case | Macro | Orbital | alpha | theta | Actual angle phi |
|---|---:|---:|---:|---:|---:|
| FeCl2 | 10 | 33 | 0.126474 | 0.197697 | 0.173732 |
| Porphine singlet | 2 | 98 | -0.565578 | 0.464440 | 0.818782 |
| Porphine triplet | 2 | 98 | -0.763049 | 0.598626 | 1.193894 |
| Porphine triplet | 2 | 100 | -0.046946 | 4.767531 | 1.373492 |

The largest porphine angles are not small perturbations of normalized additive
steps. In particular, passing theta=4.77 directly to sine/cosine is not a small
retraction correction. Negative alpha can also enlarge the actual normalized
angle: the triplet orbital 98 turns by 1.19 radians despite theta=0.60.
One rejected/shortened singlet keyframe has alpha=-1.28004 and phi=1.81987;
the normalization hemisphere can change without reaching a zero vector.

These data disprove the assumption of uniformly S-tangent raw lifts and small
per-orbital angles on these OEO trajectories. Nonzero alpha is **not** an extra
independent physical degree of freedom or, by itself, a quotient-rank bug.
The chart deliberately retains a non-spherical raw representative for full OEO;
its normalization eliminates the radial component to first order, but it changes
the finite path. Raw angles can also include removable inter-orbital/core gauge
motion. The frame-component cap does not bound these raw AO-sphere lengths.

Do not directly replace the update with exp(v) at the current numerical step.
The next relevant distinction is raw representative motion versus physical
subspace motion after core and active-gauge projection, and the relation between
the current horizontal lift and the step units. Any exponential test must account
for the changed second derivative when alpha is nonzero and retain the selected
CI transport pullback; no new empirical angular clipping is justified here.

Resource totals below include scalar/geometry logging, so are not performance
improvement claims:

| Case | Macros | HVP | Last energy (Eh) | Wall (s) | Peak RSS (GiB) |
|---|---:|---:|---:|---:|---:|
| FeCl2, converged | 17 | 322 | -2181.624464859829 | 10.02 | 0.481701 |
| Porphine singlet, short cap | 2 | 19 | -2244.7225280821 | 175.42 | 39.517925 |
| Porphine triplet, short cap | 2 | 26 | -2244.7596897644 | 167.20 | 39.533016 |

FeCl2 energy differs by +7.969e-9 Eh from completed native run 1254179. Its
17 versus 26 macros is not credited to a new optimizer: there is no optimizer
change. The first four energy rows agree to approximately 8e-9 Eh; the full
trajectory is not bitwise identical across nodes/log modes. The porphine cap
exits 2 intentionally and must not be labeled convergence.

Artifacts: `step-geometry-build-1254268`, `step-geometry-1254271`, build/run
wrappers, and `summarize_step_geometry.py`. The script exports each orbital,
macro, micro, role, stage and raw values with source-line references, plus
accepted-only CSV and aggregate CSV/JSON. Original logs remain available;
compressed logs are retained as repository evidence.

## Raw representatives versus physical active spaces (2026-10-14)

This audit uses the current finite-CI-transport implementation, not the earlier
pre-transport native trajectory. Native job 1254394 records only two macros,
with unchanged CIAH, HVP, thresholds and preconditioner. Its extra code writes
raw orbital snapshots before canonicalization; it adds no objective or HVP.
Build 1254377 passed orbital geometry, CIAH and arbitrary-selected-space CI
transport tests. The deliberate two-macro native exit is 2, not convergence.

Geometry jobs 1254391, 1254434, 1254460, 1254479 and 1254511 run on Hanhai22 compute nodes.
They use saved native/PySCF orbitals and the AO overlap only: no Hamiltonian,
many-body sigma, energy or gradient evaluation. The last job contains the full
analysis, including the local-sphere counterfactual and script/snapshot hashes.
Saved PySCF AH macro orbitals and its completed low-energy checkpoint
are reused; neither PySCF nor an old native version is rerun.

### Gauge-independent comparison

First remove the inactive span from the active columns and form an
S-orthonormal frame for that projected span. Principal angles then compare
physical subspaces, independently of raw orbital norms, phases, mixing and
export canonicalization. The saved initial native/PySCF active spaces agree
within 1.3e-6 degrees. A saved exact/raw-versus-Molden calibration bounds the
active-space angle error from eight-decimal export at 1.6e-5 degrees.

| Comparison | Largest active principal angle |
|---|---:|
| Current NEO macro 1 vs PySCF macro 1 | 5.8708 degrees |
| Current NEO macro 2 vs PySCF macro 2 | 39.7865 degrees |
| Current NEO macro 1 to macro 2 | 8.8260 degrees |
| PySCF macro 1 to macro 2 | 41.4136 degrees |
| Current NEO terminal vs PySCF low-energy terminal | 89.9997 degrees |

The old first-macro difference of about 21.94 degrees preceded production CI
transport and must not be quoted as the current result. Transport substantially
improves the first macro; physical paths diverge strongly in the second.
Terminal subspaces differ, not merely their raw representatives. Different
methods need not follow identical paths or converge to the same local minimum;
this table does not alone prove an incorrect gradient or stationary energy.
At macro 2 only one active principal direction differs strongly (39.7865
degrees); the other ten angles are below 4 degrees. At the terminal point three
directions differ by about 55.22, 55.55 and 90 degrees. It would be misleading
to describe all eleven active directions as orthogonal between the endpoints.

### Exact candidate geometry and the hidden finite-path crossing

Independent differentiation of the core and projected-active span projectors
reproduces the native accepted physical tangent norm to 1.7e-15. The second
accepted candidate and exported canonicalized active space differ by only
3.4e-6 degrees, consistent with Molden precision. Thus neither a wrong logged
step norm nor export canonicalization explains the observed physical divergence.

| Macro/stage | Physical tangent norm | Largest active angle | Active-frame polar orientation | Same-tangent isometry vs raw active angle |
|---|---:|---:|---:|---:|
| 1, accepted | 0.978722 | 38.4966 degrees | +1 | 9.0207 degrees |
| 2, first keyframe | 0.466900 | 20.9519 degrees | -1 | 21.9994 degrees |
| 2, full macro trial | 0.430968 | 30.5726 degrees | -1 | 35.0260 degrees |
| 2, half-step accepted | 0.215484 | 8.8260 degrees | +1 | 2.8532 degrees |

The isometry here is a geometry-only counterfactual with exactly the same
physical first tangent. Its skew-generator reconstruction error is below
7.2e-15. It is not an optimizer implementation, an electronic-energy test, or
permission to assume closure of an arbitrary selected VB space. Nor does it
uniformly improve a projector Taylor remainder: macro 1's remainder is larger
under that isometry, whereas macro 2's full-trial remainder is 0.151901 instead
of 0.833101 under the actual raw update.

A local-sphere exponential constructed from each normalized raw tangent was
also tested geometrically, without changing NEO or evaluating its energy. It
removes the orientation reversal in these candidates, but changes the second
keyframe's true active angle from 20.9519 to 49.9890 degrees; the full macro
trial becomes 33.1662 degrees. Thus it is not equivalent to either the raw
update or the collective same-tangent isometry. This reinforces the need to
derive the complete orbital/CI pullback before trying a new production update.

At the second keyframe, raw orbital 98 has alpha=-1.280040, theta=1.100975 and
normalized raw-ray angle 1.819871 radians (104.27 degrees). Its negative
longitudinal coefficient is real, but the resulting pathology is not just
individual-orbital scaling. Along the exact additive segment from macro origin
to this keyframe, one active principal angle reaches 90 degrees at fraction
0.781756, then the endpoint returns to an unsigned 20.9519 degrees. At that
crossing the active Gram condition is about 175.8: the active orbitals remain
linearly independent. For the full macro trial, the corresponding crossing is
at fraction 0.884797, with active Gram condition about 156.2. The path crosses
the old-active graph chart, not a singular physical active basis.

To see the mechanism, let the origin projected active columns be T0=Q0 R0,
with Q0^T S Q0=I. At a candidate let T be its active columns after removing
its current inactive span. Define W=T R0^-1, L=Q0^T S W and E=W-Q0 L. Then
Q0^T S E=0 and, whenever L is invertible,

$$
\operatorname{span}(T)
=\operatorname{span}\!\left(Q_0+E L^{-1}\right).
$$

The first projector derivative removes the active-parallel contribution, but
the finite update still depends on L^-1. In the second keyframe the first-order
active-parallel Frobenius norm is 1.350605 while the physical external tangent
norm is only 0.449683. The parallel measure includes removable column scaling;
it is not a norm of additional physical degrees of freedom. Positive scaling
does not alter a span, but its combination with an external increment changes
the finite graph denominator and can reverse its sign. The endpoint L has
smallest singular value 0.273844 and negative determinant. In the full macro
trial its smallest singular value
is 0.128225. The half-step candidate retains positive orientation.

This is why a modest physical tangent norm need not mean a modest, well-behaved
finite subspace path under the Euclidean raw lift. The complete-OEO chart uses
an Euclidean complement of gauge generators (`orbitals/charts/chart.cpp`),
whereas the physical metric uses S-projected subspace changes. These are valid
quotient coordinates at first order but do not bound the discarded parallel
representative component at a finite step.

### Coupled model failure and the PySCF difference

At the second keyframe, the actual coupled gradient norm is 0.985829 versus
model prediction 0.063051. The structure components are respectively 0.979898
and 0.053263 in the native coordinate norm. The keyframe is rejected. The full
macro trial predicts a 3.99623 mEh decrease but actually increases the relaxed
energy by 6.20311 mEh; the half step is accepted with rho=0.930389. This is a
measured failure of the finite coupled model on the candidate, not evidence
that the Newton residual tolerance should be tightened.

The native finite CI transport uses the Cholesky representation R and
representation cross-Grams R_p^T R_0 (`optimization/coupled/ci_transport.cpp`),
not the physical overlap of the two orthogonal active frames. Such a transport
is a valid coordinate retraction but is not an invariant physical-state map
under arbitrary active-frame rotations/reflections. For example an orthogonal
active rotation leaves the Gram and R unchanged, although keeping the same
physical multiconfigurational state requires an exterior-power CI coefficient
transformation. A pure active rotation is excluded from the complete-CAS
reduced chart; that example alone is not an on-chart bug. However the measured
admissible finite update above generates an actual frame-orientation reversal,
which Gram-only transport does not encode explicitly. The CI tangent is then
required to compensate a strongly nonlinear representative change.

The installed PySCF source differs concretely: `mc1step.update_rotate_matrix`
uses expmat(unpack(dx)), and `newton_casscf.extract_rotation` composes that
unitary update while adding and normalizing CI increments. Redundant
active-active rotations are not independent AH variables. On the saved PySCF
macro 1-to-2 path, the active frame polar overlap has positive orientation and
only 0.120734 degrees internal polar angle, against the native rejected
candidate's reflection. Both methods freeze an origin Hessian across keyframes;
that shared feature is not itself a missing-HVP bug. Default orbital-AH PySCF
and native coupled AH are not identical subproblem trajectories.

An orientation reversal alone is not a physical error: corresponding CI
coefficient transformations can absorb it. The issue here is its coexistence
with a large finite-path excursion, Gram-only coefficient extension and a
measured structure-gradient prediction failure. Do not equate an orbital phase
change by itself with the final energy discrepancy.

The concrete weakness now isolated is the uncontrolled finite representative
motion of the raw lift/update and its interaction with the finite CI extension.
It explains an actual second-macro model rejection and where the physical path
starts diverging. It does **not yet prove** this is the sole cause of the final
0.344725 mEh gap, or that any chosen exponential update will reach PySCF's basin.
No macro-2 derivative audit or electronic-energy test of the same-tangent
counterfactual has been made in this round. An individual-orbital sphere update
alone does not remove active-parallel mixing or supply CI frame covariance.
The next root-targeted change needs a consistent horizontal finite lift and CI
frame transport, retaining projected treatment for nonclosed selected spaces,
rather than a new clipping constant, residual threshold, or endpoint test.

The diagnostic-only two-macro run uses 19 HVP, reaches
-2244.722528082080 Eh, and has wall 222.96 s / peak RSS 39.5207 GiB. This includes
geometry logs/snapshot I/O and is not a performance benchmark or a convergence
claim. Full raw snapshots, scalar log, executable hash, test evidence and the
geometry-only script/JSON are retained in `step-geometry-1254394`,
`step-geometry-build-1254377`, and `representative-analysis-1254511`.

## Corrective S-horizontal quotient lift (2026-10-14)

The correction changes the accepted-point *linear lift*, not the orbital
normalization formula, nonlinear CI extension, HVP policy or tolerances. Let
G_p contain the existing support-admissible gauge generators on orbital p's
allowed coefficient slots and let S_p be the corresponding principal AO
overlap submatrix. The forbidden-support nullspace is constructed exactly as
before; only the complement is changed:

$$
G_p^T S_p U_p=0,
\qquad
\delta b_p=L_p z_p.
$$

Here U_p spans the nullspace and L_p is its existing normalized raw lift,
with optional coordinate whitening only for the older TN/L-BFGS charts.
Since S_p is positive definite, a vector v in both span(G_p) and its
S_p-orthogonal complement satisfies v^T S_p v=0, hence v=0. The new complement
therefore removes the same gauge rank and retains the same physical dimension.
Active--active generators still enter **only** the complete-CAS OEO chart;
arbitrary selected structures retain their independent active rays. Frozen
coefficient tails and unequal HAO supports keep their exact admissible gauge.

At a canonical full-OEO accepted point, active columns are S-orthogonal to
the core. The new lift makes core increments S-orthogonal to the old core,
and complete-CAS active increments S-orthogonal to the old occupied span.
It therefore removes the raw radial/active-parallel first jet responsible
for the measured folding. With a fixed core, if T0^T S X=0, then

$$
T_0^T S(T_0+tX)=T_0^T S T_0
$$

for every finite t. This fixes the active graph orientation for that fixed-core
path. Simultaneous core motion can still cause second-order changes; this is
not a claim of global no-reflection or guaranteed basin identity.

Changing U_p alone is insufficient. An Euclidean left inverse would no longer
annihilate the gauge. Vector recovery and covector lifting must be dual:

$$
D_p=L_p^T S_p L_p,\qquad
P_p=D_p^{-1}L_p^T S_p,\qquad
z_p=P_p\delta b_p,\qquad
g_{z,p}=L_p^T g_{b,p},\qquad
\widetilde g_{b,p}=P_p^T g_{z,p}.
$$

Consequently P_p L_p=I, P_p G_p=0, L_p^T P_p^T=I and G_p^T P_p^T=0.
This is necessary for recycled directions, transported secants and the
analytic full-OEO physical dual norm. The chart's vector projection, packed
gradient projection, gradient expansion and rotation inverse are all changed
together; no caller keeps interpreting the same inverse as an Euclidean Gram.

No per-orbital cubic inverse is introduced. For the unwhitened representative,
write n=||b||_S, c=b/n, M=U^T S U, r=U^T S c, v=M^-1 r and
delta=1-r^T v. With s=1 for the product-sphere representative and 0 otherwise,

$$
L=n(U-scr^T),\qquad
D^{-1}=\frac{1}{n^2}\left(M^{-1}+s\frac{vv^T}{\delta}\right).
$$

The existing shared M^-1 supplies this action. For a whitened B satisfying
B^T S B=I and u=B^T S c, the corresponding inverse raw AO Gram is I/n^2
for the sphere representative and (I-uu^T)/n^2 otherwise. These formulas also
cover frozen-tail charts where scaling is not admissible. The obsolete
Euclidean rank-two inverse is deleted, rather than retained as a fallback.

For the analytic full-AO preconditioner, let E=Q_ext^T S L/n be its existing
frame map. Since span(L)=span(Q_ext), its inverse is

$$
E^{-1}=n^2 D^{-1} E^T.
$$

The implementation reuses E and D^-1, and retains its existing residual
certificate; no extra AO multiplication or factorization is required.

All trial orbital coefficients are still b+Lz. L is immutable within a macro,
so the raw-orbital gradient pullback is L^T g and its Hessian pullback is
L^T H L. The finite CI extension recomputes its actual Cholesky direction using
the same L, and its gradient/curvature pullbacks remain in that coordinate
system. No derivative of an across-accepted-point basis is added to a trial,
and no exponential-retraction second derivative is silently omitted.

The opt-in keyframe derivative audit now checks
the first keyframe of **each** audited macro, so macro 2 is no longer excluded.
It is still not an optimizer acceptance criterion.

### Qualification of the implemented correction

Build 1254753 passes twelve tests: sparse-orbital quotient, arbitrary selected
CI transport, CI transport, representation overlap, overlap pullback, OEO
normalization pullback, CIAH, coupled structure, normalized orbital curvature,
analytic Hessian diagonal, orbital step geometry, and shifted-metric
preconditioner. Added assertions explicitly check S-horizontal lifts, gauge
annihilation by both vector projection and dual covector lift, unchanged
physical dimensions/supports, large fixed-core OEO finite graph orientation,
and pulled Hessian finite differences. The binary SHA-256 is
`871bb8a0f5ba78664c471f99c47537bb0e28a9b95e30c3bea3acf3ec46c49126`.

Molecular array 1254764 runs four native regressions to convergence, and the
singlet only to two macros with derivative instrumentation. All use Hanhai22,
32 bound cores, the original guesses, exact libcint, Davidson, analytic
Hessian-diagonal preconditioning and unchanged stopping thresholds. Existing
1254179 results supply the before column; no old version is rerun.

| Case | Macros before / after | HVP before / after | Wall seconds before / after | Peak RSS MiB before / after | After energy (Eh) | After-minus-before (Eh) |
|---|---:|---:|---:|---:|---:|---:|
| F2-OEO | 3 / 3 | 15 / 13 | 0.30 / 0.25 | 18.51 / 18.33 | -198.689800055684 | -5.13e-9 |
| FeCl2 | 26 / 15 | 455 / 252 | 11.69 / 6.89 | 495.82 / 518.13 | -2181.624464860406 | +7.39e-9 |
| MnF2 | 11 / 10 | 157 / 159 | 3.87 / 3.46 | 499.16 / 518.29 | -1348.893353289841 | +3.06e-9 |
| 241 | 3 / 3 | 20 / 20 | 2.01 / 2.07 | 810.52 / 843.84 | -230.720590533850 | -2.98e-12 |

These are single runs on different Hanhai22 nodes, not repeated timing
qualification. The current round targets correctness, not performance.
It does not claim universal speed or memory improvement.

Geometry-only job 1254782 analyzes the new exact snapshots with the same
projector/isometry formulas. At the first keyframe of macro 2:

| Quantity | Previous raw lift | Corrected S-horizontal lift |
|---|---:|---:|
| Max active abs(alpha) | 1.280040 | 1.49e-16 |
| Active-parallel tangent Frobenius norm | 1.350605 | 3.76e-16 |
| Active frame polar orientation | -1 | +1 |
| Actual active space versus same-tangent isometry | 21.9994 degrees | 1.1032 degrees |
| Actual structure-gradient norm | 0.979898 | 0.113134 |
| Actual / predicted coupled gradient norm | 0.985829 / 0.063051 | 0.221672 / 0.029301 |

The corrected accepted macro-2 frame has positive orientation and internal
polar angle 0.2405 degrees; its physical tangent norm agrees independently
with the native record to machine precision. Thus the previously diagnosed
raw folding mechanism is removed on the tested candidates, rather than
merely hidden by smaller step limits.

The actual first keyframe of each of the two macros also receives an opt-in
full transported HVP/gradient audit. These are derivatives of the current
finite joint map, not repeated cross-program fixed-point energy comparisons.

| Macro | Perturbation | Full HVP error norm | Finite joint energy-derivative error (Eh) |
|---|---:|---:|---:|
| 1 | 1e-3 | 1.73e-8 | 9.97e-9 |
| 1 | 3e-4 | 6.72e-9 | 1.12e-7 |
| 2 | 1e-3 | 1.66e-7 | 4.10e-8 |
| 2 | 3e-4 | 1.63e-8 | 1.07e-7 |

This round does **not** close the final singlet energy gap. The corrected
macro-2 keyframe still fails the model-mismatch test (0.229213). Its relaxed
full macro decreases energy by only 0.356948 mEh against a predicted
9.785949 mEh, rho=0.036476. Two-macro final energy is -2244.719262254443 Eh,
above the previous two-macro -2244.722528082080 Eh. The corrected and PySCF
macro-2 active spaces still differ by 41.5051 degrees. Fixing representative
drift is therefore confirmed, but neither improved early energy descent nor
identity with the final PySCF basin has been demonstrated.

The derivative-instrumented two-macro run uses 20 HVP, wall 219.22 seconds,
RSS 40469.00 MiB; its intentional exit is 2. It is not a native speed result.
Job 1254797 runs the corrected singlet to convergence from the **original**
input, with the derivative audit and raw snapshot output disabled. It retains
scalar/angle logs to detect any later recurrence; there is no new optimizer
threshold or extra endpoint criterion. Its eventual terminal energy, total
macro/HVP counts and resource totals must be reported before calling the
singlet energy issue resolved.

Evidence: `horizontal-build-1254753`, `horizontal-regression-1254764`,
`representative-analysis-1254782`, the build/run wrappers, and
`summarize_horizontal.py`. The summary contains both machine-readable CSV
and JSON, including every derivative error rather than only the best one.

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
