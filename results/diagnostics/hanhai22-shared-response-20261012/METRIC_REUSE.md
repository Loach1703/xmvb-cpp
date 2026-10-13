# Certified fixed-point metric reuse

## Regression and historical comparison

The historical reference is `results/history/paired_runs.csv`, revision
`c2776c4`, method `neo`. Use `241_vb`, not the different `241_VBSCF` guess.
All historical results are reused; no historical binary is rerun.
Hanhai25 versus Hanhai22 wall times remain contextual, not matched hardware
speedups. Historical FeCl2 reached a different, higher-energy endpoint and is
not a correct-endpoint performance reference.

Current complete structure dual norms introduced a costly Riesz solve in each
AH residual and keyframe check. In suite 1249844, QOCGUM performed 82 solves,
32,255 metric actions, and spent 41.31 s inside these solves (52.13 s total).
CERRAS performed 104 solves and 43,654 metric actions, spending 399.10 s
inside these solves (432.73 s total). These are nested stage measurements,
not sums of independent CPU times.

Scalar tracing also requested extra norms. Native validation therefore unsets
`XMVB_NEO_TRACE` and `XMVB_NEO_CANDIDATE_AUDIT`. One new-binary QOCGUM scalar
capture is separate from the native timing suite.

## Fixed-operator reuse

Let G be the unchanged positive definite structure tangent metric at one
accepted point. The physical dual norm requires the certified solution of

$$
Gx=b,\qquad \|b\|_{G^{-1}}=\sqrt{b^T x}.
$$

Previously each solve started from zero. Retain primal vectors and their
already computed metric images in a workspace belonging to one CIAH macro.
Twice-applied metric Gram--Schmidt constructs Q with Q^T G Q = I. The Galerkin
initial solution is

$$
x_0=QQ^Tb.
$$

The original positive preconditioner P is balanced against this space:

$$
P_Q=QQ^T+(I-QQ^TG)P(I-GQQ^T).
$$

This expression is symmetric positive definite: its quadratic form is the
sum of a squared coarse-space component and a positive quadratic form of
the projected component. Applications use the stored GQ and add no metric
action. In particular, b=GQ a has a known inverse image; its certification
needs one true metric action, not another CG solve.

The workspace also observes each AH direction and its existing metric image.
These pairs are available from the AH workspace without another HVP or metric
evaluation. Certified Riesz solutions provide the second family of vectors.
At most twice the existing AH direction budget is retained, independent of
unique-spin-string pair counts. Storage is O(d m), not a full structure
metric, and no per-thread copy is made. The workspace is destroyed at the end
of the macro; data are not transferred between different accepted points.

Numerical dependence controls only admission to the reuse space. It does not
change any optimization or certification tolerance. Every returned solution
still explicitly satisfies the original true-residual requirement:

$$
\|b-Gx\|_2\leq\sqrt{\epsilon_{\rm machine}}\,\|b\|_2.
$$

Complete-space factorized inverses already solve in one CG step and retain
their existing path; the reuse space accelerates selected-space Riesz solves.

## Validation ledger

- Build/unit job 1249883: initial-solution reuse; five tests passed.
- Suite 1249888: initial-solution-only prototype, original inputs and thresholds.
- Build/unit job 1249905: balanced SPD reuse; five tests passed.
- Suite 1249909: balanced reuse, 15 native systems and one new-binary scalar capture.
  QOCGUM: 3 macros / 37 HVP, 38.38 s native versus 52.13 s saved prior run;
  the separate capture used 22,548 metric actions and 29.67 s Riesz time.
  This remains much slower than the historical 6.33 s and is NOT full recovery.
- Build/unit jobs 1249932 and 1249935: also observe existing AH metric pairs;
  the final frozen build 1249935 passed all five tests.
- Final native suite 1249938: all 15 inputs submitted concurrently. A separate
  QOCGUM scalar capture reports 82 norm requests, 17,812 metric actions,
  24.21 s Riesz time, and 16 requests with no CG iterations. Native QOCGUM
  takes 33.58 s with the same 3 macros / 37 HVP and final energy within
  1e-12 Eh of suite 1249844. This is still 5.30 times historical wall time.
- Superseded prototype LOFLEA tasks 1249888_12 and 1249909_12 were cancelled;
  only the final variant's LOFLEA task is retained. Raw partial records remain.

### Final-suite snapshot (program end-to-end seconds)

| Case | Historical | Saved pre-change | New | New macros / HVP |
|---|---:|---:|---:|---:|
| F2 | 0.135 | 0.315 | 0.305 | 3 / 23 |
| 241 | 0.916 | 2.016 | 2.012 | 3 / 21 |
| 240 | 10.681 | 4.625 | 4.621 | 3 / 29 |
| MnF2 | 9.860 | 4.423 | 4.242 | 10 / 137 |
| FeCl2 | 5.285 (different endpoint) | 8.194 | 8.015 | 13 / 251 |
| 7975 | 2.059 | 3.058 | 3.050 | 3 / 25 |
| 7963 | 2.342 | 3.546 | 3.514 | 3 / 25 |
| YAMSAI | 2.091 | 1.992 | 1.910 | 3 / 26 |
| NEGPIA | 13.416 | 21.782 | 16.859 | 4 / 49 |
| QOCGUM | 6.330 | 52.133 | 33.577 | 3 / 37 |
| 67970 | 3.165 | 14.614 | 11.068 | 4 / 38 |
| CERRAS | 28.306 | 432.731 | 218.633 | 5 / 41 |
| LOFLEA | 57.215 | 1281.550 | running | - |
| 10698 | 12.352 | 14.706 | 15.573 | 4 / 46 |
| F2-OEO | 0.143 | failed | 0.199 | 3 / 15 |

The native final-suite FeCl2 energy is -2181.624464859729 Eh (within 2.5e-8
Eh of the saved strict PySCF endpoint), and peak RSS is 516.9 MiB. The full
factorized inverse path is retained. A changed macro count is not evidence
that the remaining performance issue has been fixed.

QOCGUM peak RSS is 1945.5 MiB versus the saved 1987.4 MiB; bounded reuse has
not produced a memory-growth regression in this measurement. Input initial
energies match the historical CSV within its printed precision. Remaining
regressions (including 10698's increased wall time) are explicitly retained,
not filtered from the comparison. The newly converged F2-OEO run is not a
proof that the separate zero-macro-step logic is generally repaired.

CERRAS finished with energy -397.093435003825 Eh, within 1e-12 Eh of the
saved pre-change endpoint. Peak RSS is 2017.9 MiB versus 1938.7 MiB;
the observed RSS increase must not be hidden. Its 5 macros / 41 HVP are
unchanged. Total time drops from 432.73 to 218.63 s, but remains 7.72 times
the historical 28.31 s. LOFLEA's final variant is still running at this snapshot.

Use `compare_history.py RUN_ROOT --previous SAVED_RUN_ROOT` to regenerate CSV
directly from raw outputs. The script preserves incomplete/failed runs and
compares the same program-level timer rather than GNU-time launch overhead.

The regression test covers linear combinations of previously solved RHS,
new RHS under the same certification contract, zero RHS, bounded storage,
and reuse of an already computed AH metric image.

Remaining issues are explicit: historical regressions are not all resolved,
and F2-OEO's zero-macro-step failure is separate from metric solve cost.
