# Exact selected-space metric solves and support reduction

## Scope

This change preserves the selected VB structures, orbital constraints, HVP,
NEO strategy, and convergence tolerances. It does not mark a truncated space
as complete. The baseline is the saved data in `240_subspaces.csv` and jobs
1259837/1259844; no historical executable is rerun.

## Restricted inverse

Let W_P contain the original selected spin-string expansion columns. Complete
these columns to the full target-spin sector as W = [W_P, W_Q], without
altering any entry or ordering of W_P. The completion is a metric-only object:
it neither constructs a new StructureAction nor copies electronic integrals.

For the exterior orbital transform T, define

$$
M=W^T T^T T W,\qquad N=M^{-1},\qquad G=W^T W.
$$

The full target-spin sector is invariant under the spin-independent orbital
transform and its transpose. Therefore

$$
N=G^{-1}W^T T^{-1}T^{-T}W G^{-1}.
$$

G is block diagonal by spatial occupation pattern. Its small block factors
and the sparse expansion W are shared by all solves. T is applied through
elementary exterior transforms, without forming an exterior matrix.

Partition N in the selected/complement ordering. The exact restricted inverse is

$$
M_{PP}^{-1}=N_{PP}-N_{PQ}N_{QQ}^{-1}N_{QP}.
$$

In particular, it is not N_PP. To solve M_PP x = b, compute

$$
\begin{pmatrix}y_P\\y_Q\end{pmatrix}
=N\begin{pmatrix}b\\0\end{pmatrix},\qquad
N_{QQ}z=y_Q,\qquad x=y_P-N_{PQ}z.
$$

Only the omitted coordinates are eliminated; none is admitted to the physical
variational wavefunction. This identity does not require the selected space
to be invariant under orbital transformations.

For an occupation block with d doubly occupied and m singly occupied orbitals,
the highest-weight spin dimension is

$$
\binom{m}{n_\beta-d}-\binom{m}{n_\beta-d-1}.
$$

Missing columns are obtained from the nullspace of S_+ together with the
transpose of the selected columns in that occupation block. Existing columns
remain unchanged; completion columns need not use the input generator's VB
pairing convention. Tests include unequal alpha/beta counts and column
rescaling/sign changes.

### Direct versus implicit correction

The first implementation built N_PQ and factored N_QQ at each changed orbital
chart. If q is the complement dimension and p the selected dimension, this
costs q full-inverse images plus O(q^3) factorization and O(pq+q^2) storage.
Job 1259988 confirmed accuracy, but setup dominated: ION(0:3) spent 0.477 s
of 0.523 s transport solve time constructing the correction. ION(0:2), whose
complement has 630 rather than 70 columns, regressed to 8.37 s SCF. The explicit
complement-column implementation was consequently removed, not retained as
another production branch.

The lean implementation uses P = N_PP as an SPD preconditioner, not as an
answer. With A = M_PP, B = N_PQ, and D = N_QQ, the same identity gives

$$
P=A^{-1}+BD^{-1}B^T,
\qquad A^{1/2}PA^{1/2}=I+A^{1/2}BD^{-1}B^TA^{1/2}.
$$

Only at most q eigenvalues differ from unity, and they are all at least one.
Ordinary SPD preconditioned CG solves A x = b; it does not apply ordinary CG
to the nonsymmetric product PA. In exact arithmetic it terminates after at
most min(p,q)+1 steps (and at most p by the ordinary dimension bound). This
is an algebraic bound, not a prediction of floating-point iteration count.

Each preconditioner application evaluates only N(r,0), retaining its selected
components. There are no explicit Q columns, N_PQ array, or N_QQ factorization.
The low-rank discrepancy is corrected implicitly while solving the original
selected equation. Every solution satisfies the same true selected residual
target sqrt(machine epsilon). No nested tolerance or relaxation is introduced.

Per-chart setup now only prepares the exterior transforms, shared between
retraction and gradient solves at an exactly identical chart. Memory is the
sparse completed spin expansion, occupation-block factors, and bounded carrier
workspace. A solve costs the measured number of selected metric actions plus
principal full-inverse actions; its total cost, not iteration count alone,
determines performance.

The optional completion retains a 64 MiB budget including bounded workspace.
It is not constructed if its complement is at least as large as the selected
space or its workspace does not fit. General selected-space metric CG remains
the solver for such cases, with its transformations now prepared once per
chart. This selection does not restrict admission to orthogonal-CI H/S action.
Invalid factorization or failed residual certification is an error, not an
accepted approximate answer.

## Exact support reduction

ProductPaths stores maximal generators under componentwise dominance of the
sorted occupied indices in both spin sectors. A dominated generator contributes
no additional triangular support, including its one-generator closure, and is
removed from path planning only. The actual structure expansion coefficients
are never removed or merged.

Sorting by decreasing sum of occupied indices ensures every possible strict
dominator precedes a candidate. The candidate is compared only with retained
maximal generators. With s input seeds and f retained generators, the bound is
O(s log(s) + sf(n_alpha+n_beta)); it can remain quadratic for a fully incomparable
frontier. For 240 ION(0:3), 4830 input seeds reduce exactly to two generators.

Tests compare independently enumerated supports and labels at depths zero and
one, including unequal spin counts and empty/duplicate inputs.

## Verification record

- 1259944: prepared metric and support reduction; five regression tests passed
  and all five 240 spaces converged on anode15, CPU set 64-95.
- 1259947: restricted-inverse unit test passed, but a stale test object caused
  a class-layout mismatch in selected_ci_transport. The production runs were
  correctly blocked. Preserved rsync timestamps were older than objects built
  between the two source transfers; the affected headers were touched and
  their consumers rebuilt before further testing.
- Native runs supply wall time and RSS. Separate opt-in phase runs supply
  attribution; nested timers must not be added to total time.

## Final results

Job 1260033 passed all six regression tests and completed every space with the
lean principal-inverse PCG implementation. All five native runs and the three
phase runs used one 32-core allocation on Hanhai22 anode01, CPU set
`81-95,152-167,190`. The original anode62 and then anode15 were occupied by
other long jobs; no other users' jobs were disturbed. Old results below are
saved data, not reruns. The five-space historical batch used anode62/15/24,
so old-to-new wall comparisons are not a same-node controlled speedup claim.

| Space | Structures | Steps / HVP, old and new | Old SCF (s) | New SCF (s) | New process wall (s) | New peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| COV | 14 | 3 / 20 | 2.692325 | 2.679260 | 4.28 | 1572.56 |
| ION(0:1) | 294 | 3 / 27 | 4.022312 | 3.765359 | 5.32 | 1566.61 |
| ION(0:2) | 1134 | 4 / 32 | 5.680133 | 4.515104 | 6.03 | 1606.10 |
| ION(0:3) | 1694 | 4 / 32 | 7.588366 | 4.691147 | 6.34 | 1622.50 |
| FULL | 1764 | 4 / 30 | 3.602792 | 3.545873 | 5.15 | 1572.07 |

The saved controlled old ION(0:3)/FULL pair on anode62 had SCF times
6.033106/3.397351 s and process times 7.67/5.11 s. Thus the within-node
ION(0:3)/FULL SCF ratio shrank from 1.78 in that saved run to 1.32 in the new
same-allocation run. The baseline CSV's 9.29 s ION(0:3) process time should
not be substituted for the controlled pair without noting provenance.

| Space | Final energy (Eh) | New minus saved old (Eh) |
|---|---:|---:|
| COV | -343.239211501262 | approximately -1e-12 |
| ION(0:1) | -343.407736636063 | approximately +2e-12 |
| ION(0:2) | -343.443415326853 | approximately -2e-12 |
| ION(0:3) | -343.446256564846 | approximately -1e-12 |
| FULL | -343.446302545202 | approximately -1e-12 |

All initial energies match the saved inputs/results to printed precision.
No convergence tolerance or HVP rule changed. Peak process RSS did not
decrease materially: the problem is dominated by other allocations. The new
prepared metric preconditioner reports 151996 bytes of persistent storage for
both ION(0:2) and ION(0:3); it is shared, not copied per thread.

Separate phase runs quantify the work reduction:

| Transport metric quantity | Saved old ION(0:3), 1259844 | New ION(0:3), 1260033 | New ION(0:2), 1260033 |
|---|---:|---:|---:|
| Solves | 38 | 38 | 36 |
| Selected metric actions | 7058 | 688 | 940 |
| PCG iterations | 7020 | 650 | 904 |
| True residual certifications | 38 | 38 | 36 |
| Solve wall (s), including preparation | 2.252324498 | 0.575635797 | 0.763119928 |
| Maximum certified relative residual | 1.4863633e-8 | 1.48156862e-8 | 1.47108504e-8 |

Each new PCG iteration additionally applies the completed principal inverse;
the 10.8-fold iteration reduction is not a 10.8-fold wall speedup. New orbital
transform setup totals 0.000547/0.000549 s for ION(0:3)/(0:2), respectively.
The selected solve remains exact to the original numerical target.

`240_metric_results.csv` contains native results for the prepared-only,
explicit-complement diagnostic, and final lean implementations. It can be
regenerated from saved run folders with `tests/mcscf/summarize_240_metric.py`.
