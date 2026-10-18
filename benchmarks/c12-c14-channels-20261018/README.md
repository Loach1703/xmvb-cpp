# Native C12/C14 follow-up

The user switched performance validation from C16 to C12/C14. New C16
jobs293254/293255 were cancelled before completion; their partial timing is
not an end-to-end result. Historical executables are not rerun.

Frozen executable: qualification293253,6/6 tests in6.57s, SHA256
`141283957128192fe63fac1a82120b1cb98202876b06242bf97bad5f2a81a7f5`.
It includes compact beta target rows and fused exchange weighting.

Use original archived inputs, exact Libcint/cc-pVDZ, CIAH/Davidson,32 physical
cores, OpenMP close/cores, BLAS1. All cases run to their existing convergence
criteria and original input iteration limits. No extra SCF step cap is set.
C12 here is the established CAS(12,12) LOFLEA benchmark, not a carbon-C12
molecule. C14 is anthracene CAS(14,14),COV429.

Cases are C12 COV132,ION(0:2)47256,FULL226512 and C14 COV429. The main runs
use the current4GiB shared CI budget. A fifth run checks C12 COV with64MiB,
the archived budget. This tests resource scaling without rerunning old code.
Archived runs used Hanhai22 EPYC9654; new runs use Hanhai25 EPYC9755. Runtime
ratios across clusters or budgets cannot isolate this kernel change.

Relevant saved controls:

| Case | Job | Steps/HVP | Application s | GNU process s | GNU RSS GiB |
|---|---|---:|---:|---:|---:|
| C12 COV |1264892|4/36|17.208014|17.40|2.5721|
| C12 ION(0:2) |1265266_1|4/44|104.286|107.11|2.824|
| C12 FULL |1265524_2|4/40|100.682823|111.18|3.5782|
| C14 COV |1265059|4/46|306.172876|306.62|5.7971|

Older C14 control1264476 took1190.925493s application time,4/46,energy
-535.672671742972Eh. It must not replace the better306.17s saved control.
See `benchmarks/history/c12_c14_comparison.csv` and
`results/diagnostics/c12-amplitude-20261017/README.md` for frozen evidence.
The FULL control is the better saved1265524_2 measurement, rather than the
older1265114_2 measurement (102.727418s application/112.47s process).

Array293257_[0-4] submits these five cases concurrently on Hanhai25.
Native results and numerical comparisons will be recorded after completion.
Build time is outside the native resource timer. Frozen input SHA256 values:

| Case | SHA256 |
|---|---|
| C12 COV |092d545498678b0cf772c5c84b9eab674b694375558cd391a1b739b7bdcc35eb|
| C12 ION(0:2) |48b44d64b773d4abaf388fcf728d478a09ab76f5426f889822b1dbc8e031432f|
| C12 FULL |3f612eab8d7d74fc6e13f1ffe1c9b96864366d664c4dadc7e96dcd4e155287b9|
| C14 COV |c9600592898c7cadf087d9cb3f77854a83eddcccc51962ac3c20bb0f393520b3|

## Completed results

All five native runs exit0 with `ciah_projected_gradient_tolerance`. The
original input hashes and initial energies match their saved controls.
`current.csv` is parsed from scalar output records, CIAH traces and GNU time;
`comparison.csv` joins it to the saved historical tables by input hash.

| Case | Steps old/new | HVP old/new | Application s old/new | Process s old/new | GNU RSS GiB old/new |
|---|---:|---:|---:|---:|---:|
| C12 COV |4/4|36/33|17.208/15.227|17.40/15.42|2.572/3.423|
| C12 ION(0:2) |4/4|44/50|104.286/81.984|107.11/82.66|2.824/3.848|
| C12 FULL |4/5|40/44|100.683/106.673|111.18/109.18|3.578/4.387|
| C14 COV |4/4|46/39|306.173/219.998|306.62/220.36|5.797/7.985|

New energies are -399.106646267319,-399.397232509245,
-399.417495587054 and -535.672671743023Eh, respectively. The largest
difference from the selected controls is4.8e-11Eh. HVP/trajectory counts
change, so these are end-to-end observations, not equal-action kernel ratios.

C12 COV's additional64MiB run293257_4 converges in4 steps/33 HVPs,18.985039s
application/19.12s process,2.972900GiB RSS,energy -399.106646267317Eh.
Both budgets use the same frozen executable but different nodes. Larger
shared residency improves this observed runtime at increased RSS; the main
4GiB results are not evidence of a memory reduction or an isolated change
versus historical64MiB runs. No large numerical array is replicated per worker.

Raw native artifacts are retained remotely under
`/home/guqqgroup/taoxia/work/xmvb-ciah-subspace-20261017/benchmarks/c12-c14-channels-20261018/run-293257`.
Local diagnostic records retain traces, resource logs, affinity and hashes in
`results/diagnostics/c12-c14-channels-20261018/run-293257/`.

## Correct current-RHS replay after the FULL step regression

The five-step FULL run cannot stop at step4: its energy decrease is
9.7894656e-7Eh, above the unchanged1e-7Eh threshold. Historical FULL reached
its near-final energy during macro3 with22 HVPs; unconditional clipped-path
replay instead used11 and left a correction for macro4 and a fifth acceptance.

The driver formerly consumed an unused AH path without testing its current
right-hand side. An exact source Ritz residual does not certify that
remaining direction after the gradient changes. Use cached H/G images to
form its complete current AH residual, including outside-space components,
and apply the existing1/3 forcing fraction. Failure discards only the unused
path and reprojects in the retained span before any enrichment. No new HVP,
metric solve, true gradient, threshold or molecule-dependent policy is added.
The formula and independent augmented-matrix check are recorded in
`article/ciah_clipped_path.md`.

Qualification293280 passes8/8 tests in7.23s, including `ciah` and `keyframe`.
Frozen binary SHA256:
`ab08632d8da2d6be5952760eda828dc6bd0b96b21d29a4c6de19a12bd44c1e5e`.
Array293281_[0-3] converges all four original inputs on Hanhai25 with32 cores,
the same4GiB shared budget and unchanged outer criteria. Scalar event traces
do not enable physical-model audits. No prior executable is rerun.

| Case | Steps before/after | HVP before/after | Process s before/after | GNU RSS GiB before/after |
|---|---:|---:|---:|---:|
| C12 COV |4/4|33/33|15.42/15.12|3.423/3.345|
| C12 ION(0:2) |4/4|50/44|82.66/77.29|3.848/3.869|
| C12 FULL |5/4|44/40|109.18/104.58|4.387/4.503|
| C14 COV |4/4|39/46|220.36/222.50|7.985/7.986|

COV has zero replay resets and preserves its trajectory. The other three
cases each reject three stale remainders. FULL's per-macro HVP counts return
to8,7,22,3 and its energy is -399.417495587054Eh, unchanged from the preceding
run. ION's energy is unchanged; C14 ends at -535.672671742975Eh (4.8e-11Eh
from the preceding result). C14's HVP count increases; this correction is
not a demonstrated all-system reduction of action counts. COV/ION/FULL use
the same respective nodes as293257; C14 changes from anode060 to anode033.
All these nodes are EPYC9755, but the allocations are not exclusive repeats.

The two frozen measurements remain separate in `current.csv` and `replay.csv`;
`replay-comparison.csv` pairs their original input hashes. Raw logs and traces
are in the corresponding `run-293281/` directory alongside `run-293257/`.

## End-to-end audit: unnecessary inverses in final structure normalization

The scalar CI/HVP/keyframe records are not a complete runtime partition.
FULL's process time is104.58s, while the recorded CI solves, HVP expansions
and fixed-structure keyframes account for17.35s,9.31s and7.34s, respectively.
The norm of every individual structure is still required for the normalized
structure-coefficient report; Coulson--Chirgwin weights additionally use the
retained M*C images. Neither requires a spin-overlap inverse at reporting time.

Audit293293 enables existing phase scopes on the then-current executable
from qualification293280. Its process time is80.44s on anode014; four CIAH
macros total34.25s. This is attribution, not a matched speedup measurement.
Sampling293294 on anode019 resolves the call chain
`build_exact_structure_overlap_diagonal -> exact_spin_overlap -> resolve_matrix
-> FullPivLU inverse -> dtrsm -> BLAS workspace mutex`.
GNU gprofng2.42 reports a changed collection timer; its numeric sample shares
are not reliable wall fractions. The resolved unnecessary inverse call is
also directly established by the source, independently of that warning.

For each pair of expansion terms, a structure norm consumes only the product
of the alpha and beta overlap determinants, multiplied by their two expansion
coefficients. The inverse matrix, cofactor data and signed-log determinant
state are unused. The determinant-only entry point therefore shares the
existing factorization and LU/SVD rank decisions but omits inverse solves and
unconsumed matrix copies/state. This is not a regular-pair approximation and
does not change the singular-value threshold. Empty, signed, rank-deficient
and relative-LU/absolute-SVD rank cases are checked in `cofactor_differential`;
`orthogonal_direct_ci` checks structure norms against an independently formed
small dense overlap matrix for arbitrary expansion coefficients.

The physical gradient/HVP, all CIAH policies and outer tolerances are
unchanged. The existing structure-local overlap maps remain local and are
cleared between structures: this change adds no global pair cache, no complete
structure matrix, and no worker copies of the integrals. It does not remove
the remaining pair-of-expansion-terms work; cross-structure reuse is a separate
possible improvement. `structure_normalization` is an opt-in phase around the
actual final calculation, not an additional calculation or certificate.

`normalize.sbatch` builds/tests only the changed exact implementation and
freezes its source manifest and binary. `audit.sbatch` records native phase
timing or gprofng sampling on a frozen executable. Native performance
comparisons use the existing293281 records, without rerunning historical
executables. Unrelated local RI edits are excluded from the remote mirror.

### Completed qualification and native comparison

Qualification293303 passes10/10 tests in6.67s, including singular/cofactor,
arbitrary-structure overlap, CI transport, CIAH, and the F2 text report.
Binary SHA256:
`5a7d1d68b7e7bab5706f5f35d12b4d3666f3e8d6a227f13d98af7b5665c7a681`.
Native array293313_[0-3] converges every original input with32 physical cores,
the same4GiB shared CI workspace, scalar tracing and unchanged thresholds.

| Case | Steps/HVP before = after | Energy after (Eh) | Process s before/after | GNU RSS GiB before/after |
|---|---:|---:|---:|---:|
| C12 COV |4/33|-399.106646267317|15.12/13.30|3.345/3.454|
| C12 ION(0:2) |4/44|-399.397232509243|77.29/67.45|3.869/3.859|
| C12 FULL |4/40|-399.417495587056|104.58/47.09|4.503/4.541|
| C14 COV |4/46|-535.672671742974|222.50/212.15|7.986/7.980|

Maximum final-energy difference from293281 is2.0e-12Eh at printed precision.
Each per-macro HVP sequence is also preserved. RSS is not systematically
reduced: this change removes computation, not a dominant resident matrix.
FULL's complete process speedup is2.22x; its application timer changes from
102.094402s to45.003713s. These are single measurements on different nodes
of the same EPYC9755 cluster; small timing changes are not isolated kernel
speedups. Historical Hanhai22/64MiB records remain separate.

New native phase audit293314 measures `structure_normalization` once at
0.267026176s. Its four CIAH macros total34.008920629s, compared with
34.250658471s in audit293293; removing reporting-time inverses did not
accelerate the physical macro actions. Physical-metric solves total only
0.054133184s in this audit (0.046387482s before): repeated physical-norm
reporting is not the missing large cost here and was not given a new cache.

The remaining main costs are different: ION still needs273 accumulated CI
actions versus FULL's103, and C14 remains dominated by CI/HVP work. The
structure-local pair contractions remain exact and are not claimed to have
a lower asymptotic exponent in this round.

`normalization.csv` and `normalization-comparison.csv` are generated by
`summarize.py` from saved report tails, CIAH traces, hashes and GNU-time
records. It checks input hashes, successful convergence and energy agreement
before joining with `replay.csv`. Raw large output files remain on Hanhai25;
compact logs, report tails and frozen build metadata are mirrored in
`results/diagnostics/c12-c14-channels-20261018/run-293313/`.
