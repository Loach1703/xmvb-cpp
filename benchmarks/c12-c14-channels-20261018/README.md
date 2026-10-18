# Native C12/C14 follow-up

## Trial-chart keyframe overlap reuse

The transport point now retains shared immutable ownership of its own exact
selected Gram, when that Gram already exists. The fixed-state keyframe consumes
it instead of replaying an overlap action. The normalized overlap image is
`(M*C)*white`; the sigma/integral-adjoint sweep receives this known image and
omits the second structure-overlap projection and the primal reread used only
for that projection. Hamiltonian images, integral cotangents, generator moments,
physical structure response and all solver policies/thresholds are unchanged.
Selections without a resident Gram still compute the required overlap, then
reuse its normalized image. No Gram is added solely for this optimization.

The whole tile regression was added to qualification. Initial293355 passed
12/13 tests but exposed a boundary bug from the preceding column-plan change:
sampled spin-string membership queries were incorrectly required to provide
a complete alpha carrier. Prefix membership is now kept separate from the
complete-carrier column generator in `mcscf/ci/covalent_columns.hpp`. A sampled
membership view does not create a full universe or a product table. The
original high-orbital-count prefix tests are retained, not weakened.

Final qualification293374 passes13/13 tests in5.87s, including structure tiles,
fixed-state energies/gradients/integral adjoints with and without a supplied
Gram, known-overlap consumption, unequal spins, arbitrary selected structures,
multistate transport and immutable trial-Gram lifetime. Frozen binary SHA256:
`58b4c9b76d787aa5db031b16d4a12c3a2f301430441fe877a7b0dcbf399a740c`.
The failed qualification and successful replacement are both retained.

Native293375_0 (C12) and293376_3 (C14) run concurrently on Hanhai25 with the
same exact inputs,32 physical cores,4GiB shared CI workspace, CIAH/Davidson,
BLAS1 and original outer criteria. No archived executable is rerun.

| Case | Steps/HVP old = new | Final energy new (Eh) | Best saved / preceding / new process s | RSS preceding / new GiB |
|---|---:|---:|---:|---:|
| C12 COV132 |4/33|-399.106646267318|13.30 /17.39 /13.80|3.403 /3.429|
| C14 COV429 |4/46|-535.672671742974|212.15 /202.49 /201.80|7.947 /7.857|

Initial energies and per-macro HVP sequences match. Final differences from the
saved excellent controls are at most1.1e-12Eh. GNU time measures the complete
application process; application-report timers are also retained in the CSV.
C14's new time is4.9% below the excellent control but only0.34% below the
preceding version: this round is not a demonstrated significant C14 speedup.
C12 still does not beat13.30s.

There is concrete placement evidence relevant to C12's preceding slower run:
excellent293313_0 used CPUs96-127 (one NUMA domain); preceding293347_0 used
96-105,116-127,224-233 (two NUMA domains). New293375_0 uses160-191 (one NUMA
domain) on anode003 rather than anode031. Therefore its17.39->13.80s change
cannot be assigned solely to these source changes. This allocation evidence
does not by itself quantify a NUMA penalty. C14 remains on anode064, but its
core mask changes as well. All these nodes report EPYC9755 processors.
No binding policy, thread count or physics parameter was changed to obtain
the new results; the actual affinity files are retained.

Current-only audit293377, launched after the C12 native run, verifies that
pure-overlap product actions fall9->0 while Hamiltonian actions remain79.
Compared with saved audit293348, fixed-state keyframes now use the transported
trial Gram and no longer replay the overlap at a changed chart. Normalized
overlap reuse also removes the product-jet primal-overlap projection. C14's
remaining native major scopes are81.44s CI,80.05s physical HVP and13.92s
fixed-state keyframes: further large end-to-end gains must target the actual
CI/HVP contraction, not claim another large win from the now-removed metric
replays.

`keyframe.csv` and `keyframe-comparison.csv` compare against the excellent
normalization293313 controls. `summarize.py --case-job c14-cov=293376` records
the separate array allocations accurately. Source/build/test manifests,
native resource logs, actual affinity and phase evidence are retained under
the qualification/run/audit job directories in
`results/diagnostics/c12-c14-channels-20261018/` and the remote benchmark mirror.

## Selected metric and symbolic column reuse (current follow-up)

This round changes neither the Hamiltonian/HVP nor CIAH policies, physical
spaces, thresholds, workspace allowance or thread count. The coupled CI metric
now uses the already-resident exact accepted-point selected Gram when present.
Selections without that Gram retain their exact streamed metric action.

Complete covalent support prepares a coefficient-independent prefix DAG once
per support object. Its three edge/subtree flags use
`U_beta * (n_act+1) * (n_alpha_electrons+1)` bytes, with linear single-spin
permutation/prefix data. This is not a product-sized mask or numerical cache.
Complete subtrees emit colex intervals; a temporary one-column bitset restores
arbitrary caller index order without sorting. CI and derivative tiles use the
same row-generation interface and immutable reachable beta list. General
selected structures retain their exact maximal-generator constraints: they
are not incorrectly given the complete-covalent support certificate.
The derivation is in `article/exterior_column_schedule.md`.

Qualification293346 passes12/12 tests in1.80s. Tests independently check
prefix support for unequal spins and permuted indices, arbitrary selections,
the explicit horizontal metric, absence of metric-triggered product generation,
multistate transport derivatives and existing orthogonal-CI actions.
Frozen binary SHA256:
`d500edd47116a4af6bd3cfb8ca1dfaecbb8e8c92fd691a8253b098a1b19b7d80`.
`projection.sbatch` retains configuration, source hashes, build and test logs.

Native array293347_[0,3] runs the two unchanged inputs to convergence with
32 physical cores, exact Libcint/cc-pVDZ, CIAH/Davidson, BLAS1 and the same
4GiB shared CI budget. No archived executable is rerun. The comparison uses
the better saved normalization293313 records, not the slower intermediate
or cross-cluster controls. Times below are GNU-time complete process times;
application report timer boundaries changed during the output-layer refactor.

| Case | Steps old/new | HVP old/new | Energy new (Eh) | Process s old/new | RSS GiB old/new |
|---|---:|---:|---:|---:|---:|
| C12 COV132 |4/4|33/33|-399.106646267318|13.30/17.39|3.454/3.403|
| C14 COV429 |4/4|46/46|-535.672671742975|212.15/202.49|7.980/7.947|

Both initial energies match and final differences are1.1e-12Eh or less. Each
per-macro HVP sequence is unchanged. C14 is4.6% faster on the same anode064
EPYC9755 node. C12 uses anode031 instead of the archived anode019, both EPYC9755,
and does **not** beat the13.30s excellent record. Its CI time increases from
3.855s to5.499s despite the same79 CI actions; HVP time also increases from
3.738s to4.809s despite the same33 expansions. The present evidence does not
attribute that difference or demonstrate an all-system end-to-end speedup.
Against the immediately preceding saved output-layer COV run293330_0,
17.85s/4steps/33HVP, the new native time is17.39s, a small observed improvement.

Current-only phase audit293348 verifies the actual removed work. Compared with
saved audit293333, pure-overlap product actions decrease33->9: the coupled
accepted-point metric replays disappear, while the nine trial-chart keyframe
normalizations remain. They cannot use an origin-chart Gram. The79 Hamiltonian
actions are unchanged; their cumulative `ci_product_pattern` scope falls from
1.353187s to0.307037s (4.4x). Including the removed metric calls, that scope
falls1.865668s->0.343368s. These are diagnostic scope measurements, not an
end-to-end speedup or additive wall-time partition.

Initial qualification293341 and native293343/audit293344 measured the then-new
DAG before replacing permuted-index sorts with bitset emission; they are retained
as intermediate evidence, not substituted for final qualification293346 and
native293347. No historical baseline was rerun in either implementation check.
`projection.csv` and `projection-comparison.csv` are parsed by `summarize.py`
using `--cases c12-cov c14-cov`; raw reports and resource logs remain in the
corresponding remote `run-293347/`, with compact local diagnostic copies.

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

## Output-only structure normalization

The optimizer no longer completes exact structure norms for reporting. The
objective's reporting-only method and its call are deleted, without an alias.
The result retains only the small accepted active-orbital overlap needed by
optional analysis, not an HVP context or copied physical pair cache. It is
important to use this accepted overlap: exported OEO orbitals can be in a
different representative from the stored structure coefficients.

`output/text/normalization.*` owns the output preparation. It reuses an
already exact result diagonal (or a matrix already present) before doing any
new contraction. Otherwise it prepares topology once and calls the existing
overlap-only kernel. Normalized coefficients, Renormalized weights and the
reported mean overlap share one norm vector across all states. Energy,
optimization and orbital writers do not request it. Coulson--Chirgwin weights
still use the retained raw C and M*C and do not require structure norms.
No full H/M matrix or eigensolve is added; Davidson remains the API/input/CLI
default and is checked explicitly. Existing dense data are merely reused if
an explicit reference calculation has already produced them.

The SCF timer now excludes optional report normalization, while the command
end-to-end timer is sampled after final structure analysis, including its
output. Consequently application timers across this change have different
boundaries; GNU process time remains the comparison metric.

Qualification293329 passes11/11 tests in1.92s. Tests cover arbitrary selected
open-shell expansions, unequal alpha/beta electron counts, accepted-frame
overlap without AO/integral input, exact-diagonal/matrix reuse, coefficient
and weight invariance, Davidson defaults, and a two-state F2 report with at
most one norm contraction. An initial build293326 reached the application
but failed to find the new test target in Make's stale top-level rules; the
script now explicitly reconfigures before building. Binary SHA256:
`dfb30ac85e7896c6044921709b4f785e32c51bc43b9b3a3aebc985872313b082`.

Native293330_0/2 preserve4/33 and4/40 steps/HVPs for C12 COV/FULL. Every
printed normalized coefficient, Coulson--Chirgwin weight and Renormalized
weight agrees exactly at the printed8 decimals with saved293313 output:
132 entries per table for COV and226512 per table for FULL. No baseline
executable is rerun. `output-layer.csv` records the complete native results.
FULL takes46.23s process time,4.597GiB RSS,energy-399.417495587054Eh;
COV takes17.85s,3.372GiB,energy-399.106646267316Eh. This refactor is not a
claimed physical-kernel speedup. The COV native measurement is slower than
293313: its79 CI actions take5.59s rather than3.85s, and33 HVPs take4.92s
rather than3.74s. Nodes differ; the unchanged counts alone do not explain
that timing difference or establish its cause.

Current-only phase audits293331(FULL) and293333(COV) each record exactly one
output norm contraction,0.355720514s and0.007424226s, respectively. The FULL
audit's macro work totals33.81s; COV's totals12.29s. These scopes separate
output work from the physical solve and show that the COV timing increase
is not seconds of new structure normalization. No new threshold, regular
certificate, physical action or per-state normalization is introduced.
