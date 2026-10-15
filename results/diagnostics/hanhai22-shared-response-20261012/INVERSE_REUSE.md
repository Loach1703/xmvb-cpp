# CIAH inverse-image reuse

## Scope

The macro-point coordinate chart and metric remain fixed. This change retains
the inverse image and full linear defect of the expensive structure-metric
solve. It does not alter CIAH, HVP, coordinates, preconditioner, or convergence
thresholds. The cheap analytic orbital dual norm is unchanged.

For the signed AH reference component a, the residual is

    r = a (g + H p - omega G p).

Writing x_g = inverse(G) g and e_g = g - G x_g, the clipped model update is

    g_new = g + t H p
    x_new = (1-t) x_g + (t/a) x_r + t omega p
    e_new = (1-t) e_g + (t/a) e_r.

At a true keyframe, the mismatch inverse and defect are the differences of
the true and predicted gradient inverse images and defects. These are full
space identities, not projected AH residuals. Negative a must retain its sign.

## Implementation

- `MetricImage` holds the inverse and full linear defect.
- One cached RHS permits retrieval of the residual inverse already evaluated
  by CIAH norm certification without another metric action.
- Clipped model gradients and keyframe mismatches use the algebraic inverse
  as the initial solution. One true metric action verifies each new RHS;
  CG only refines it if the existing linear-residual tolerance is not met.
- The refinement target remains sqrt(machine epsilon) times the **new RHS**
  Euclidean norm. It is not scaled by the correction RHS or the parent RHS.
- The propagated defect is retained for diagnostics, not trusted as a substitute
  for the true residual after floating-point cancellation.
- Unclipped model norms retain the existing AH scalar identity. Their inverse
  images are propagated lazily until needed. Gradient rollback also restores
  the corresponding inverse image.

This eliminates repeated solves, not every metric application. A small
Euclidean linear defect is not a condition-independent bound on physical
dual-norm error; no stronger certification claim is made here.

Persistent additional storage consists of a bounded number of structure-sized
vectors. No structure matrix, carrier array, or per-thread large array is added.

## Verification

New unit cases exercise signed roots, clipped/unclipped steps, full-space
defects, exact RHS cache hits, and near-cancelling keyframe gradients. A
deliberately incorrect advertised defect must not bypass true-residual
refinement. An exact supplied inverse must finish with one metric action and
zero CG iterations.

Build/unit job: 1250146. New-binary regression array: 1250148 (Hanhai22,
32 bound cores). Fifteen native cases plus one independent QOCGUM scalar
trace. Old versions are not rerun. The immediate saved reference is
`norm-decisions-final-1250047`; historical reference is revision c2776c4 in
`results/history/paired_runs.csv` (Hanhai25, different cluster).

All five unit executables passed, including the new inverse-reuse tests.

QOCGUM's new scalar capture records 31 exact-RHS cache hits, 12 seeded
combinations, and 1086 CG refinement iterations across those combinations
(1.594 s). None of the 12 inherited combinations met the unchanged relative
linear-residual target without refinement. Total recorded structure-norm
work is 13959 CG iterations / 14067 metric actions / 20.298 s.

The preceding capture records 12943 iterations / 13031 actions / 18.290 s,
but **omits the old mismatch predicate's internal solves**. These totals
therefore cannot be used to claim that inverse reuse increased or decreased
total metric actions. Native QOCGUM time is 39.806 s versus saved 33.170 s;
the current result is not an end-to-end speedup. Both use 3 macros / 37 HVPs
and reach the same final energy. Timings are single runs on different nodes.

Fourteen native cases completed and converged; LOFLEA remains running at this
snapshot. Completed final energies differ from the immediately preceding run
by at most 5.41e-9 Eh (FeCl2). CERRAS retains 5 macros / 41 HVPs, with
183.605 s versus 179.871 s and 1969.7 versus 1962.0 MiB peak RSS.
FeCl2 retains 12 macros, with 249 versus 248 HVPs. F2-OEO uses 3 versus 4
macros, with 15 HVPs in both. Other completed macro/HVP counts are unchanged.
See `riesz-reuse-comparison.csv` for all native results and the historical
comparison, including slower cases. This implementation is not evidence of
universal speedup or restored historical performance.
