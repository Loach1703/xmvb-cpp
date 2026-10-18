# C16 CIAH startup diagnosis

## Scope and outcome

Tracked source: `aa222faa`, with no production changes. The unrelated local
RI edits were not transferred: the execution mirror was refreshed from a Git
archive. The original attachment is `C16_vb.xmi`; its bytes, geometry, basis,
guess and tolerances were preserved. It declares ISCF=8, so the diagnostic
explicitly passes `--optimizer-backend ciah`, as requested by the user.

This is C16H18, neutral singlet, CAS(16e,16o), STR=COV, strictly localized
HAO orbitals, exact Libcint/cc-pVDZ, Davidson, 330 Cartesian AOs and 1,430
VB structures. The execution confirms 12,870 strings in each spin sector.
Davidson remains in the selected VB structure space, not in product space.

The bounded test reproduces the apparent hang **before the initial CI solve
has returned**, not in an AH macroiteration. The principal sampled work is
the right exterior transform used to generate orthogonal-CI input panels.
No convergence failure or incorrect energy has been established by this test;
there is no initial/final energy or completed orbital-HVP count to report.

## Reproduction and resource evidence

Hanhai25 Slurm **292840**, anode017, AMD EPYC 9755, 32 physical cores,
affinity 224-255, OMP close/cores, BLAS one thread. Release C++23/GCC14.3,
`-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`. Source, binary and input hashes,
build flags and environment are retained in the raw records.

| Native bounded run | Observed |
|---|---:|
| Wall / s | 120.24 |
| User CPU / s | 3094.39 |
| System CPU / s | 641.11 |
| CPU utilization, including waits | 3106% |
| Peak RSS, GNU time / KiB | 7510076 |
| Peak RSS / GiB | 7.1622 |
| Minor page faults | 43454829 |
| Major page faults / swaps | 2 / 0 |
| Exit | 124, intentional diagnostic timeout |
| Completed outer steps | 0; output ends at the iteration header |

The job's successful batch exit only means the diagnostic capture completed.
It does not mean the molecular calculation converged. The 120-second limit
belongs to the diagnostic wrapper, not a changed SCF threshold or ITMAX.
Busy/waiting OpenMP workers are included in the CPU count; this is not a
measurement of 31 useful computational cores.

A separate 75-second gprofng 2.42 clock capture, including startup, resolves:

1. `ExteriorOrbitalTransform::transform_right` task body, the leading sampled
   application function;
2. `gomp_team_barrier_wait_end`, substantial waiting;
3. `ExteriorOrbitalTransform::apply_supported_left`;
4. support generation, sorting and projection, smaller in this capture.

The selected sample CPU shares are approximately 49.6%, 36.6% and 7.1% for
the first three entries. They are **not stage wall-time percentages** and
must not be converted into a speedup prediction. The collector follows the
`stdbuf` exec and records about 74 seconds for the actual program. It reports
no errors or timer-change warnings, but the program's child experiment is
marked not closed because the bounded collection was terminated. Use this
as call-path evidence, not a complete-run timing or allocation profile.

Trace confirms direct-CI support construction and the spin metric
preconditioner (one 1430-dimensional block, 16,359,200 factor bytes), but no
first Davidson Ritz-iteration trace or CIAH keyframe is reached. The observed
path is initial structure CI -> product-tile sigma input reader -> panel
generation -> right exterior transform. Thus neither repeated AH expansion
nor orbital-structure response is the evidenced startup bottleneck.

## Symbolic scale and the residency cliff

`support_counts.py` counts the complete-covalent prefix paths used by
`ProductPaths`. This is combinatorial planning, not an electronic-structure
run. For half-filled singlets, each initial product has complementary alpha
and beta occupations. After an upper-triangular exterior transform the
symbolic support obeys the prefix condition

$$
N_\alpha(k)+N_\beta(k)\ge k.
$$

The code's depth-one closure uses the corresponding bound with $k-1$.
These are planned traversal supports, not numbers of independent physical
degrees of freedom or guarantees that every floating-point amplitude is
nonzero. The selected variational dimension remains the VB structure count.

| CAS | Structures | Strings per spin | Symbolic support | Depth-one closure | Input payload / MiB | Accepted payload / MiB |
|---|---:|---:|---:|---:|---:|---:|
| 12,12 | 132 | 924 | 226512 | 608751 | 2.59 | 4.32 |
| 14,14 | 429 | 3432 | 2760615 | 7770620 | 31.59 | 52.65 |
| 16,16 | 1430 | 12870 | 34763300 | 101508836 | 397.83 | 663.06 |

Input payload is one double and one int row index per supported product;
accepted payload is two doubles (Y and its residual image) sharing one int
index, excluding small column offsets/ownership metadata. These are payload
estimates, **not measured total RSS**. The C14 row is a dimension comparison,
not a rerun or a matched comparison of different molecules' wall times.

Three concrete execution issues follow from the current implementation:

1. `orthogonal_tiles.cpp` quarters a 256 MiB workspace and tests full-product
   residency. For C16 a dense vector is 1263.71 MiB, so RHS batches are one
   and the ordered resident representation cannot be used. A compact input
   also exceeds its 64 MiB share. H/S source panels therefore cannot all be
   retained during the sigma neighbor sweep.
2. `ColumnPanels::get` rebuilds a missing input panel and consumes/projects
   a partial output panel when evicted. `ExteriorOrbitalTransform::row_block`
   separately has a 64 MiB LRU and regenerates a missed row block by applying
   the exterior adjoint to unit rows. Sigma's neighboring beta columns can
   revisit many panels. Thus capacity is bounded, but transform/projection
   work is not bounded by one visit per panel. This is the repeated path
   supported by the sampling capture; exact miss counts are not yet recorded.
3. `product_jet.cpp` retains accepted Y/residual panels only if the entire
   accepted state fits and no output eviction occurred. C14's payload is
   near the present budget; C16 cannot qualify. Later HVPs would lose the
   anchored reuse of accepted images and replay work. This third issue is
   a source-derived consequence, **not a measured C16 HVP result**: the
   startup test never reaches an orbital HVP.

`apply_supported_left` also scans the full shear link list for every beta
column and tests allowed source/target rows inside that loop. Restricting
stored support does not by itself remove these link visits. Sorting/packing
and regenerated panels add work beyond the required Slater-Condon channels.

## What this establishes, and the next implementation target

The current exact arbitrary-selected-structure projection is not proved
wrong by this test. However, a correct projected action and a low-memory
tile implementation do **not** imply cost proportional to the number of
selected structures. Product support grows much faster here, and independent
LRU layers repeatedly rebuild its representation. The current broad-support
execution strategy is inadequate for this C16 startup.

The next task is the contraction schedule, not changing CIAH tolerances:
measure row-block/input/output rebuild multiplicities, then organize the
shared exterior maps and sigma neighbor traversal as one bounded contraction
plan. Physical overlap must remain exact, and source/projection reuse must
be accounted across consumers. It is not valid to remove the exterior maps
by assuming closure of the selected COV space. Nor should simply increasing
the fixed 64 MiB constant be reported as an algorithmic solution. An exact
low-memory algorithm has a time/storage tradeoff that the new plan must make
explicit; one cannot promise zero replay at arbitrary size with fixed memory.

No production code, threshold, physical model, optimizer policy or supplied
orbital guess was changed. No old executable was run. No long-running C16
job remains after this bounded diagnostic.

## Artifacts

- Local raw logs: `results/diagnostics/c16-startup-20261018/run-292840/`.
- Remote logs, executable and complete gprofng experiment:
  `/home/guqqgroup/taoxia/work/xmvb-ciah-subspace-20261017/benchmarks/c16-startup-20261018/run-292840/`.
- `run.sbatch`: bounded native run plus one clock-sampling capture.
- `support_counts.py`: reproducible symbolic counts and payload estimates.

The superseded submission 292833 was cancelled before running a molecule
when the user clarified CIAH. Its interrupted compile left three zero-byte
linalg object files; 292834 and 292839 failed to link, not to calculate C16.
292840 removed/rebuilt only those generated linalg products and linked
successfully. The unsuccessful attempts are not runtime evidence.
