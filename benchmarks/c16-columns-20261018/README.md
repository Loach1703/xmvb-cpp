# C16 exact-demand contraction repair

## Scope and result

This changes the structure H/S product-tile contraction schedule, not CIAH's
physical model, optimizer policy, preconditioner or accuracy thresholds.
Davidson still acts in the selected VB structure space. No old executable
was rerun. The unrelated local RI changes were not transferred to the execution
mirror and are not part of this repair.

The final 32-core short run reaches the first two-vector CI Ritz iteration,
but **does not finish the initial CI solve in 120 seconds**. No converged CI
energy, molecular total energy, orbital HVP count or convergence speedup is
claimed. The exact-demand schedule reduces bounded-run peak RSS by 31.0%.
The remaining real column revisits are still expensive.

## Mathematical changes

See `article/exterior_column_schedule.md` for the derivation.

1. Generate only beta exterior rows actually requested by sigma, rather than
   their enclosing contiguous panels. Separate batched generation from
   single-column residency. Consume cache hits before admitting misses.
2. Generate an exterior row by descending occupied-row Laplace recursion,
   enumerating only triangularly allowed combinations. Two single-spin arrays
   per worker suffice; no inverse, dense compound matrix or U-squared table
   is needed. Colex predecessor ranks use shared prefix/suffix sums.
3. Intersect allowed-alpha and shear-source bitsets instead of scanning every
   rejected link. Share the same graph traversal across RHS lanes while
   retaining shear and edge order within each RHS.
4. Retain only actually written output columns. Move evicted partial columns
   into a bounded projection batch. Linearity permits exact projection of
   partial images; unwritten enclosing-panel columns need no work.
5. Borrow a budget-fitting group of output owners and update their disjoint
   beta/RHS/alpha ranges in one task group. This avoids one OpenMP join per
   column without copying integrals or CI arrays.
6. Accumulate overlap once per input column per RHS block, even if the column
   is regenerated. RHS width depends on sigma workspace, not on whether the
   entire spin-product grid fits in memory.

The existing cache/workspace budgets were not increased. A 64 MiB per-worker
scratch target limits row-generation concurrency; one worker is the unavoidable
minimum if two single-spin arrays exceed that target. At C16/32 workers this
scratch is about 6.28 MiB. Pending output projections have the existing bounded
batch extent, not an unbounded list of evictions.

## Matched short-run evidence

Original attachment: C16H18, singlet, CAS(16,16), COV, HAO, Libcint/cc-pVDZ,
330 Cartesian AOs, 1,430 structures, 12,870 strings per spin. Input bytes are
unchanged. Its ISCF=8 is overridden with `--optimizer-backend ciah` explicitly.

Both measurements use Hanhai25 **anode017**, AMD EPYC 9755, 32 physical cores,
affinity 224-255, OMP close/cores, BLAS one thread; GCC14.3/C++23 Release
`-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`. Times exclude compilation. The
120-second timeout is a diagnostic wrapper, not an SCF iteration or accuracy
limit. GNU time includes startup and any work before interruption.

| Native bounded measurement | Stored baseline 292840 | New 292879 |
|---|---:|---:|
| Wall / s | 120.24 | 120.14 |
| Peak RSS / GiB | 7.1622 | 4.9440 |
| User CPU / s | 3094.39 | 3713.30 |
| System CPU / s | 641.11 | 65.70 |
| Minor page faults | 43454829 | 39613464 |
| First CI solve completed | No | No |
| Completed orbital steps | 0 | 0 |
| Wrapper exit | 124 | 124 |

**Logging qualification:** the old native run did not enable CI iteration
traces. Absence of its native Ritz messages is therefore not evidence that
the first Ritz iteration takes more than 120 seconds. Its separate 75-second
traced sampling run did not return a Ritz iteration. The new native run did
enable the existing scalar trace and returns iteration 1 (basis 2). These
different progress windows do not establish a full-action or CI speed ratio.
CPU totals include OpenMP waits, not only useful arithmetic.

The new first Ritz active-CI energy is -20.772251659301592 Eh with relative
residual 0.0015795880413146854. This is an **unconverged active-space Ritz
value**, not the molecular total energy. The separately instrumented run
returns -20.772251659301553 Eh and residual 0.0015795880413148012. Their
agreement is a consistency observation, not an independent C16 energy oracle.

## Complete first sigma sweep of the instrumented run

The final phase run completes the initial two-RHS H/S action before starting
the next correction action. The sigma executor reports:

| Stage | Accumulated wall / s |
|---|---:|
| Prepare channels | 14.725686521 |
| Read input columns | 37.717152107 |
| Compute | 7.126228890 |
| Write / partial projection | 36.028737642 |
| Sum | 95.597805160 |

This is an instrumented sigma-sweep stage sum, **not** the full initial CI
solve, nor a native H/S wall time: final projection flushes outside sigma and
other CI setup are excluded. The sweep has 10,069 batches, 50,756,849 targets,
795,230 requested source columns and a 67,105,144-byte peak sigma payload.
Input/output eviction counts are 231,694/233,364. Relative to 12,870 unique
beta columns this is still roughly eighteen evictions per column. Eliminating
unrequested panel members does not eliminate these genuine graph revisits.

During development, the same first 256 batches (371,635 targets) exposed
write-stage costs of 70.948950872 s (292868), 1.798166121 s after exact-demand
output projection (292875), and 0.037528056 s after shared writer task groups
(292879). These prefixes defer different amounts of final projection work;
**they must not be advertised as complete-action speedup factors**. The final
whole-sweep breakdown above is the useful remaining-cost evidence.

## Numerical qualification and artifacts

Final numerical code passes 5/5 selected CTests on a compute node: ci_transport,
selected_ci_transport, orthogonal_direct_ci, structure_tiles and
coupled_structure (2.48 s with OMP=2 in 292879). Tests include independent minor
determinants for every particle sector at n=3,6,8, reversed string order,
noncontiguous/repeated requested rows, empty requests, execution inside an
existing OpenMP team, and a forced tiny input-column cache. Existing exact
H/S oracles, arbitrary selections, unequal spins, exchange and derivative
checks remain active. Final include/comment/indent cleanup passes the same
5/5 tests in build-only Slurm job **292889** (2.45 s, anode017); the frozen
source manifest is verified and the executable relinked. C16 is not rerun.

- Baseline report: `benchmarks/c16-startup-20261018/README.md`.
- Machine-readable native comparison: `results.csv`.
- Local raw logs: `results/diagnostics/c16-columns-20261018/run-292879/`.
- Remote frozen executable, logs and manifests:
  `/home/guqqgroup/taoxia/work/xmvb-ciah-subspace-20261017/benchmarks/c16-columns-20261018/run-292879/`.
- `run.sbatch` builds/tests then captures native and phase runs; setting
  `C16_BUILD_ONLY=1` performs compilation, tests and artifact freezing only.
- Binary SHA256: `3d3d3a11b62bef4bd2d828366f4f77966e8a939112c3f5b77a109104710ca87e`.
- Input SHA256: `044a228603ede8cf877b72889f71e65e5d7d7a3d4e1ec950fd75f609ddd0d448`.
- Basis SHA256: `d53e7e02d55fb99cdf823d50cf13fec375e0cbe1220ff61816aee7ab38e38ec9`.

The final numerical job ends normally as a diagnostic batch, with both
molecular wrappers intentionally returning 124. No long C16 run is left
running. The separate accepted-point jet/HVP residency is **not solved or
benchmarked** by this startup change. The next target is scheduling true
source-column reuse against projection output lifetimes; increasing the
budget alone is not an algorithmic solution.
