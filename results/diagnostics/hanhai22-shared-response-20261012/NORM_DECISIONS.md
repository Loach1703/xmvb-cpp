# Progress estimates and certified CIAH decisions

Theory: [ciah_norm_decisions.md](../../../article/ciah_norm_decisions.md).

No CIAH, HVP, coordinate, outer-convergence, or mismatch threshold was changed.
No full physical metric or thread-local large matrix was added. The norm
comparison uses a variational lower bound on the complete covector; unresolved
comparisons retain the original true-residual Riesz solve. No estimated global
spectral lower bound is used to claim a small residual.

Changes:

1. Defer physical norm evaluation before the existing AH start cycle unless
   the eigenvalue-change prerequisite makes convergence possible.
2. Certify direction dependence using the complete residual, not its AH
   projection. Include the uncaptured residual in the trial inverse image.
3. Cache each candidate's complete norm and derive its model-gradient norm
   from the same solve, retaining the signed cross term.
4. Certify a keyframe mismatch threshold rather than computing a discarded
   exact norm when a variational lower bound already rejects it.
5. Replace diagnostic-only norm solves with explicitly labelled estimates.
6. Distinguish model agreement from real-energy acceptance. Save already
   evaluated descending keyframes independently of model agreement; the
   existing macro energy check still decides whether to accept them.

The last item fixes a specific failure observed in F2-OEO. In the scalar
capture 1250044, macro 2's first microstep predicted a gradient norm
0.000164579182, while the true keyframe norm was 0.000905510372. The unchanged
mismatch criterion rejected continued use of the model, although the true
keyframe energy had decreased to -198.689797786773 Eh. Rolling back the only
microstep discarded this candidate and produced `ciah_zero_macro_step`.
The repair preserves the energy candidate, not the stale model.

## Evidence

- Cluster: Hanhai22, one task, 32 bound cores, BLAS one thread, exact libcint,
  Davidson, analytic Hessian diagonal, original inputs and stopping thresholds.
- First build/unit tests: 1250016, all five passed.
- Preliminary native array and separate QOCGUM scalar capture: 1250017.
- New-binary F2-OEO failure capture: 1250044. No historical binary rerun.
- Final build/unit tests including the energy-checkpoint repair: 1250046,
  all five passed.
- Final parallel array: 1250047, 15 native systems and one separate QOCGUM
  scalar capture. Native times, not scalar-capture times, establish performance.

All 15 final native systems have now converged. CERRAS completed in
179.871131 s (5 macros / 41 HVPs, -397.093435003823 Eh, 1961.97 MiB RSS);
LOFLEA completed in 448.997780 s (4 macros / 34 HVPs,
-399.106646007264 Eh, 2817.65 MiB RSS). The previously failing F2-OEO converges in
4 macros / 15 HVPs, to -198.689800055623 Eh, in 0.223930 s with 18.296875 MiB
GNU-time peak RSS. This confirms the energy-candidate repair without a
threshold change. Its prior saved converged endpoint differs by 6.6e-11 Eh.

Selected final native results (application end-to-end time; one run each):

| System | Historical time (s) | Previous saved time (s) | Final time (s) | Macros / HVP | Energy (Eh) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| 241 | 0.915561 | 2.011902 | 1.983050 | 3 / 21 | -230.720590392829 | 828.84 |
| MnF2 | 9.859770 | 4.242068 | 4.984673 | 10 / 137 | -1348.893353223758 | 530.29 |
| FeCl2 | 5.285014* | 8.014698 | 7.923583 | 12 / 248 | -2181.624464872126 | 499.46 |
| QOCGUM | 6.329732 | 33.577303 | 33.170255 | 3 / 37 | -480.094766681464 | 2011.32 |
| 67970 | 3.165093 | 11.067634 | 10.366835 | 4 / 38 | -442.247231417072 | 1956.89 |
| 10698 | 12.352040 | 15.573286 | 13.759796 | 4 / 46 | -422.611824623576 | 4982.71 |
| CERRAS | 28.305941 | 218.633093 | 179.871131 | 5 / 41 | -397.093435003823 | 1961.97 |
| LOFLEA | 57.214518 | 559.903089 | 448.997780 | 4 / 34 | -399.106646007264 | 2817.65 |

* Historical FeCl2 is not the same final-energy endpoint. Differences between
native nodes and clusters remain; these measurements do not establish a
universal improvement. Final MnF2 and YAMSAI are slower than the immediately
preceding saved runs. All completed final rows are retained in
`norm-decisions-final-1250047/historical-comparison.csv`.

The preliminary CSV deliberately retains the failed F2-OEO and pending LOFLEA
rows; it is not the final-repair result. Most small-system iteration/action
counts are unchanged. Preliminary CERRAS is 184.961694 s versus the saved
previous 218.633093 s, both 5 macros and 41 orbital HVPs; energy differs by
approximately 2e-12 Eh. This does **not** close the gap to its historical
28.305941 s result.

Compare final results with the saved, already-run prior version and historical
CSV, without executing either old version:

```bash
python3 compare_history.py norm-decisions-final-1250047 \
  --previous metric-reuse-1249938
```

Historical `results/history/paired_runs.csv` is revision c2776c4 on Hanhai25;
current results are Hanhai22. FeCl2's historical endpoint is higher by about
6.825 mEh and is not an equal-correctness performance baseline. Raw output,
GNU-time RSS, input/binary hashes, node information, and action traces are
retained in each run directory.

## Remaining limit

Re-derivation (local manuscript note):
`article/ciah_norm_rederivation.md`. The next concrete optimization is to retain
complete Riesz images and propagate them through clipped gradients and keyframe
differences. At a fixed macro metric these are linear combinations, not
independent inverse problems. Their full residual defects must be propagated
as well; a small difference RHS may require refinement under the original
linear-residual criterion. No production implementation of this follow-up was
made during the derivation review.

Terminology correction: certifying a small Euclidean linear-system residual
does not give a condition-independent error certificate for the physical dual
norm. A known preimage in the orthogonal spin-product carrier supplies a separate
upper bound without estimating the global minimum metric eigenvalue; for an
arbitrary selected structure space this bound need not be exact or tight.

Large-residual rejection can be cheap, but proving a small residual still
requires the complete solve unless a reliable full-space upper bound is
available. Lower projected norms or estimated Ritz minima cannot supply that
certificate. These changes remove premature and reporting-only work; they do
not claim to eliminate the large selected-space metric cost.
