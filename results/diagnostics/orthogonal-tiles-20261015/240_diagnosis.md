# 240: truncated ionic space versus FULL

## Controlled measurement

Current implementation only; no historical executable was rerun. Jobs 1259837
and 1259844 each held one 32-core allocation on Hanhai22 anode62, CPU set
136-167, for all cases. Native runs and phase-instrumented runs were separate.
Job 1259844 added only opt-in CI-transport solve counters, not an algorithm,
tolerance, basis, initial-orbital, or memory-budget change. Its two CI transport
regression tests passed. Both spaces converge in four outer iterations, with
32 and 30 orbital HVPs respectively.

| Native run, job 1259844 | ION(0:3) | FULL |
|---|---:|---:|
| Structures | 1694 | 1764 |
| Initial energy (Eh) | -343.075601975739 | -343.075645280883 |
| Final energy (Eh) | -343.446256564845 | -343.446302545202 |
| Input preparation (s) | 1.498100 | 1.563773 |
| SCF (s) | 6.033106 | 3.397351 |
| Application wall (s) | 7.538410 | 4.968763 |
| Process wall, GNU time (s) | 7.67 | 5.11 |
| Process peak RSS (KiB) | 1598372 | 1680460 |

The first same-node run, 1259837, gave 7.96/5.15 s process wall and
6.355/3.509 s SCF. Thus the earlier different-node comparison exaggerated the
gap, but did not cause it. Final energies agree with the saved earlier results
to about 1e-11 Eh; the two different physical structure spaces need not have
identical energies.

## Attribution

The following is job 1259844 with `XMVB_PHASE_TIMING=1`. Nested phases must
not be added to their parents or to total time.

| Phase / counter | ION(0:3) | FULL |
|---|---:|---:|
| SCF wall (s) | 6.441321 | 3.413957 |
| CI transport metric solves | 38 | 32 |
| CI transport solve wall (s) | 2.252324498 | 0.022564693 |
| Metric CG iterations | 7020 | not used |
| Metric CG actions | 7058 | not used |
| True-residual certifications | 38 | not measured by this counter |
| Time inside CG metric actions (s) | 2.223439252 | not used |
| Keyframes: calls / seconds | 13 / 2.430357 | 10 / 0.676560 |
| Orbital HVP: calls / seconds | 32 / 0.559408 | 30 / 0.523008 |
| Sigma: calls / seconds | 455 / 0.425541 | 228 / 0.131679 |
| Left transforms: calls / seconds | 15638 / 0.418871 | 1846 / 0.023606 |
| Right transforms: calls / seconds | 14587 / 0.730902 | 1259 / 0.059721 |

The transport-solve excess alone is about 74% of the instrumented SCF excess.
This does not establish that all of it is removable. The selected solver uses
184.7 CG iterations per solve on average, at most 256. It performs exactly one
true-residual certification per solve; repeated certifications are not the
problem here. The maximum reported relative residual is 1.4863633e-8, with the
unchanged sqrt(machine-epsilon) stopping rule. Zero residual values in the
factorized-mode log are unused sentinels, not measured residual certificates.

`CiTransport::solve` calls `RepresentationOverlap::apply(upper, upper, x)` for
each CG action. That function rebuilds four exterior-transform objects each
time, expands to the carrier, applies two forward and two adjoint transforms,
then contracts to structure coefficients. For 7058 actions this accounts for
28232 transform-object constructions and 14116 left plus 14116 right transform
calls. The latter explains almost all of the observed right-transform count.
The operator and upper matrix are fixed during each solve: rebuilding those
objects is unnecessary. In this singlet, alpha and beta strings are identical,
so the four immutable transformations also represent the same object.

FULL instead uses `StructureMetricInverse::apply_at_transform`, an exact
factorization in a complete spin sector. The selected path has no analogous
efficient restricted solve; it starts fresh diagonally preconditioned CG for
every RHS. This is distinct from the H/S tile kernel and explains why tuning
only beta batches cannot fix most of this end-to-end gap.

Neither perf nor gprofng is available on this compute-node PATH (saved probe
files). Attribution is from opt-in boundary timers and exact action counters,
not sample percentages. Transform-construction time is included in the metric
action total, not separately measured.

## Two exact algebraic opportunities

### 1. Restricted metric inverse

Let the complete VB Gram matrix be partitioned into retained structures P
and omitted structures Q. ION(0:3) omits 70 fully paired spatial occupations;
the corresponding singlet structures have no unpaired-spin multiplicity.
For a verified matching FULL structure basis, write

$$
M=\begin{pmatrix}M_{PP}&M_{PQ}\\M_{QP}&M_{QQ}\end{pmatrix},\qquad N=M^{-1}.
$$

The required inverse is **not** the PP block of the full inverse. Instead,

$$
M_{PP}^{-1}=N_{PP}-N_{PQ}N_{QQ}^{-1}N_{QP}.
$$

For RHS b, form y=N(b,0), solve N_QQ z=y_Q, and return y_P-N_PQ z.
The correction sets the omitted coefficient components to zero; it does not
allow omitted structures to relax and does not assume selected-space closure.
This is a 70-dimensional exact complement correction in this case, not an
approximate physical expansion. It is a candidate, not yet implemented or
timed. Basis ordering/signs must be verified, and building 70 inverse columns
at every changed upper matrix must be compared with the measured CG work.
Storage can remain bounded by applying these columns in blocks. For arbitrary
small selected spaces, the complement may be large, so this is not a universal
constant-cost solution.

The immediate exact reuse opportunity is simpler: prepare the metric action's
immutable transformations once per upper matrix, not once per CG iteration.
Reuse between retraction and gradient solves is valid only at an identical
upper matrix. No tolerance change is required by either route.

### 2. Canonical generators of triangular support

Actual production diagnostics report 70 alpha strings, 70 beta strings and
4830 distinct original products for ION(0:3), versus 4900 for FULL. At half
filling the truncated seeds are exactly pairs a != b. The 70 excluded original
products have four doubly occupied orbitals.

Under the existing upper-triangular exterior transform, a seed generates its
componentwise downward closure in sorted occupied-orbital indices. In zero-based
indexing let T=(4,5,6,7), the unique maximum, and S=(3,5,6,7), its unique coatom.
Every string except T lies below S. Consequently the two seeds (T,S) and (S,T)
generate the same closure as all 4830 seeds. The exact transformed support has
4899 products (only (T,T) is absent); its one-generator closure has all 4900.

The standalone diagnostic uses the actual `ProductSupport` implementation for
both seed sets and verifies identical emitted sets at depth 0 and depth 1.
Hanhai22 job 1259846:

| Traversal | Products | Original bitset word operations | Two-generator word operations |
|---|---:|---:|---:|
| Depth 0 | 4899 | 753160 | 9910 |
| Depth 1 | 4900 | 2240480 | 29480 |

That is exactly 76 times fewer bitset word operations, **not** a measured
76-fold wall-time speedup. Current `ProductPaths` only sorts/deduplicates seeds;
it does not remove dominated generators. The same seed tests recur in input
transforms, H/S back-projections, output masks and sigma path construction.

In general, a seed dominated in both spin sectors contributes no new depth-0
or depth-1 support and can be removed exactly. An efficient frontier constructor
is needed; a quadratic all-seed dominance scan could replace one bottleneck
with another. This simplifies an exact path representation, not the physical
VB subspace or its coefficients. The original 1694 structure vectors and the
final structure projection must remain unchanged.

## Priority

1. Eliminate invariant transform reconstruction inside transport metric solves.
2. Evaluate the exact restricted inverse/complement formula against current CG,
   including setup, multiple RHS, and bounded memory; preserve the original residual.
3. Canonicalize triangular support generators and recognize exact full halos;
   do not enlarge the physical selected structure space.

This investigation changed diagnostics only. No optimizer rule, physical
coordinate, integral, convergence tolerance, or production support algorithm
has been altered.
