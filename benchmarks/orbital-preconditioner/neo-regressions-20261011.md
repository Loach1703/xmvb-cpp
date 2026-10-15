# NEO regression audit and optimization

## Evidence and scope

The reference is the current-version panel from job `281497`, saved in
`results/history/neo_validation_20261011.csv`. The older comparison remains
the `c2776c4` NEO rows in `results/history/paired_runs.csv`.

All scientific executions and builds ran in Hanhai25 Slurm compute jobs.
Each panel used the same 22 saved historical inputs, 32 compact-affinity
OpenMP cores, one BLAS thread, exact Libcint integrals and Davidson.
The projected-gradient and adjacent-energy tolerances remain 1e-3 and
1e-7 Eh. Additional production stationarity and symmetry checks are retained.
No molecule-specific thresholds, dense Hessian or MPI replication were added.

Jobs and immutable executable hashes are recorded with each native run:

- `281537`: opt-in stage timers on the reference executable for 241,
  LOFLEA and NEGPIA. Timers are nested and must not be added together.
- `281541`: candidate-policy changes and five unit tests.
- `281543`: candidate-policy-only ablation on all 22 systems.
- `281565`: final build and seven unit tests, plus MnF2 native/trace checks.
- `281566`: final native 22-system panel.

Comparisons are single native runs, not same-node paired repetitions; modest
timing and RSS differences can include load and allocator variability. These
measurements do not establish a universal speedup for every possible input.

## Findings

### AO preparation was paying for a global sort

Reference 241 spent 0.851 s building AO integral input, versus 0.129 s in the
historical run, despite similar integral-provider times. Reference 7975 had
the same pattern: 1.440 s versus historical 0.211 s.

The current graph builder sorted all primary integral-edge references for
locality. Preserve exactly that storage order, but produce it by scanning CSR
rows instead. Primary edges are injective, so temporarily encode their column
indices as negative markers, count primary references per row, and write each
row's references into its disjoint output range. Restore the column indices
before returning the graph. This preserves values, duplicate input records,
pair orientation and ascending edge order.

The reference-ordering work changes from

$$
O(N_{\mathrm{ERI}}\log N_{\mathrm{ERI}})
\quad\text{to}\quad
O(N_{\mathrm{edges}}+N_{\mathrm{AO-pair}}).
$$

For symmetric CSR, $N_{\mathrm{edges}}\leq 2N_{\mathrm{ERI}}$.
Only one shared $O(N_{\mathrm{AO-pair}})$ prefix buffer is added; existing
integral-index buffers are overwritten in place. Integral values are not
copied and there is no integral-sized per-thread workspace. The compact
large-operator representation is unchanged.

241's measured stages isolate this effect:

| Variant | AO input s | Preparation s | SCF s | Total s |
| --- | ---: | ---: | ---: | ---: |
| Reference 281497 | 0.851 | 1.131 | 0.871 | 2.006 |
| Candidate policy only 281543 | 0.839 | 1.122 | 0.807 | 1.932 |
| Candidate policy + linear reference ordering 281566 | 0.157 | 0.431 | 0.800 | 1.236 |

### Exact structure certification preceded candidate readiness

In coupled NEO, a projected candidate satisfying the raw physical KKT target
was sent to the exact residual-response solver before the code checked whether
it was ready for an actual nonlinear trial. This charged an exact response
even when the caller already allowed an uncertified energy/gradient keyframe.

The reference stage audit measured 241's six structure solves at 0.181 s,
with 0.046 s of structure-adjoint work. LOFLEA's six structure solves took
37.632 s, with 0.961 s of structure-adjoint work. These are measured stage
times, not sample-derived percentages.

Move the existing candidate-readiness test before this certification. It
returns `TrialReady`, not `Converged`. Actual relaxed energy and gradient
still validate the candidate, and the caller still disables early trials
near stationarity. Rejected early candidates resume the retained workspace.

For a genuinely certified step, the existing residual-response check remains:

$$
r_o=g+Ap+B^Tq+\lambda Mp,\qquad r_s=Bp+Cq,
$$

$$
\delta q=-C^\dagger r_s,\qquad
\widetilde r_o=r_o+B^T\delta q,\qquad
\widetilde r_s=r_s+C\delta q.
$$

Here $A$, $B$, $C$ and $M$ denote supplied orbital, coupling, structure and
metric operators, not explicitly assembled full matrices. The corrected
physical residual is checked to the same requested accuracy. Nothing about
the expensive certification is weakened for a step labeled `Converged`.

### Core-only candidate policy was not connected

The production driver passed an early-trial permission to coupled NEO but
discarded it in core-only NEO. Connect the same coordinate-residual candidate
heuristic to the core solver. Its default remains a certified solve; explicit
curvature diagnostics and hard cases cannot take the early branch.

This is a consistency fix, not an explanation of all NEGPIA regression:
the ablation reduces NEGPIA from 170 to 168 HVPs but leaves 17 accepted points.
The reference timer puts 12.175 s in its 170 HVPs. Its nonlinear trajectory,
not structure-response certification, remains a separate unsolved issue.

## Final panel

22/22 runs report convergence. All 22 final energies agree with the
pre-change current panel within 4.365e-9 Eh. None increases accepted-point or
orbital-HVP counts. All 22 measured wall times are lower than that panel;
the median reference/current ratio is 1.447. This is an observed cross-run
comparison, not a controlled same-node estimate of a 44.7% speedup.

Each cell below is **pre-change current / optimized current**.

| System | Accepted points | HVP | Structure actions | Total s | RSS MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| 101687105 | 5/5 | 12/12 | 157/96 | 1.756/1.137 | 960.0/981.4 |
| 104448 | 7/7 | 41/41 | 275/170 | 10.073/5.952 | 3500.4/3572.8 |
| 14069533 | 6/6 | 38/38 | 221/79 | 16.938/10.811 | 5894.8/5957.7 |
| 240 | 8/8 | 13/12 | 0/0 | 3.650/2.398 | 1728.5/1704.4 |
| 241_VBSCF | 8/8 | 27/27 | 151/87 | 2.634/1.822 | 933.7/918.6 |
| 241 | 6/6 | 10/10 | 119/79 | 2.006/1.236 | 930.5/958.7 |
| 67970 | 7/7 | 13/13 | 0/0 | 4.820/3.328 | 2009.6/1953.0 |
| 7000 | 7/7 | 42/42 | 277/171 | 8.292/5.504 | 3307.8/3354.5 |
| 7963 | 8/8 | 14/14 | 158/138 | 3.960/2.220 | 1679.4/1727.6 |
| 7975 | 5/5 | 11/11 | 100/38 | 2.961/1.665 | 1615.3/1573.1 |
| F2 | 5/5 | 11/11 | 20/5 | 0.181/0.149 | 60.0/62.0 |
| F2_OEO | 3/3 | 19/19 | 17/7 | 0.196/0.152 | 62.0/60.0 |
| F2_SA2 | 5/5 | 10/10 | 25/5 | 0.152/0.117 | 62.0/60.0 |
| FeCl2 | 3/3 | 287/287 | 1008/1008 | 5.944/5.877 | 416.8/425.4 |
| MnF2 | 14/14 | 180/180 | 487/437 | 5.017/4.027 | 483.9/480.9 |
| CERRAS | 6/6 | 30/30 | 619/331 | 36.556/30.868 | 2026.6/2035.4 |
| CONYAH | 7/7 | 13/13 | 0/0 | 3.125/1.926 | 1648.9/1613.4 |
| LOFLEA | 6/6 | 25/25 | 393/202 | 102.028/84.596 | 2822.6/2808.7 |
| NEGPIA | 17/17 | 170/168 | 0/0 | 22.798/17.167 | 3770.2/3837.4 |
| QOCGUM | 7/7 | 15/15 | 0/0 | 8.997/7.515 | 2100.2/2129.2 |
| YAMSAI | 5/5 | 12/12 | 156/60 | 1.684/1.017 | 861.0/853.5 |
| YEPKAF | 7/7 | 25/25 | 0/0 | 12.505/9.014 | 2346.3/2286.5 |

RSS is native GNU-time process peak memory. Small increases are present;
there is no claim of uniformly lower measured RSS. The code-level new
storage is one shared AO-pair prefix array, not large per-thread replication.

Full initial/final energies, hashes, historic NEO and historic block-LBFGS
records, and pre-change current references are retained in
`results/history/neo_optimized_20261011.csv`.

## Remaining regressions

Compared with `c2776c4`, observed faster samples increase from 7/22 to 13/22.
Nine remain slower. Historical and current final stationarity criteria are
not identical, so the old timings are not an equal-accuracy oracle; final
energy discrepancies from that revision remain explicitly in the CSV.

The next distinct investigations are:

1. NEGPIA: 17 points versus six historical accepted points; connecting the
   early-trial heuristic is insufficient. Inspect nonlinear model/step
   quality and curvature differences on identical accepted-point orbitals.
2. FeCl2: 287 HVPs and 1008 structure actions remain unchanged; the long
   physically certified second step dominates, not the removed sort.
3. LOFLEA/CERRAS: certification still incurs 202/331 structure actions.
   Reuse or improve the response solve without dropping its accuracy check.

No universal convergence improvement or quadratic-rate result is claimed.
This round removes a broadly shared preparation overhead and avoids premature
certification, while preserving tested energies and stopping accuracy.

## Verification

Job `281565` passed `test_neo`, `test_neo_response_solver`,
`test_neo_globalization`, `test_sparse_orbital_quotient`,
`test_shifted_metric_preconditioner`, `test_active_two_electron_sparse`, and
`test_ri_ao_h1e_hvp`. Added tests verify early trials cannot claim a physical
certificate, strict continuation retains certification, and linear integral
ordering preserves reversed pairs, duplicate records, diagonal/empty rows and
all column indices. No old sorting or uncertified-solve fallback was retained.
