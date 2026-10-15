# Same-node phase and memory audit — 2026-10-10

## Scope and controls

Baseline: a013272bd1a763679faab94b7370261f4302cd36 plus diagnostic timers.
Updated: 254554d plus the same timers and a projector payload counter.
No numerical algorithms or thresholds changed during this audit.

Hanhai25 jobs: parallel builds/tests 280379; four paired native cases and
separate heaptrack captures 280381; allocator ablation 280385; saved-capture
analysis 280387. All completed successfully. Both builds use the same
dependency environment and Release `-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`,
with Eigen BLAS enabled. Each molecular pair runs baseline then updated on the
same exclusive node and 32-core allocation, with one BLAS thread, exact
libcint integrals, Davidson, gradient tolerance 1e-3 and energy tolerance
1e-7 hartree. Independent pairs run as parallel array tasks. These are single
paired measurements, not statistically established speedups; updated runs are
always second, so warm-up/order effects are not eliminated.

Raw input/binary hashes, node/CPU/environment, outputs, native resource logs,
phase counters, and heap reports are retained locally under
`results/studies/orbital-phase-audit-20261010/` and remotely under
`/home/guqqgroup/taoxia/xmvb-runs/orbital-phase-audit-20261010/`.

## Native results

GNU-time process wall time and RSS (not batch elapsed or profiled RSS):

| System | Method | Steps, old → new | Wall s, old → new | RSS MiB, old → new |
|---|---|---:|---:|---:|
| 241 (HAO), anode022 | block-LBFGS | 9 → 9 | 1.20 → 1.12 | 962.1 → 969.6 |
| 241 (HAO), anode050 | NEO | 6 → 6 | 1.59 → 1.51 | 967.7 → 972.1 |
| FeCl2 (OEO), anode055 | block-LBFGS | 256 → 255 | 67.82 → 62.93 | 428.0 → 428.1 |
| FeCl2 (OEO), anode075 | NEO | 3 → 3 | 11.66 → 11.56 | 449.2 → 448.4 |

All converged. Energy differences: zero at printed precision for 241/block-
LBFGS; approximately 1e-12 Eh for both NEO cases; 6.442e-9 Eh for FeCl2/block-
LBFGS. The earlier conspicuous slowdown is not reproduced in these controls.
The earlier cross-node, concurrent-job comparison cannot isolate an actual
implementation regression from execution variability.

## Why the isolated projector speedup is not the total speedup

The 7.15-times result was a synthetic 469-AO/89-inactive/11-active **projection
only** benchmark. It excluded chart construction, curvature eigendecompositions,
integrals, CI, HVPs, and every other optimizer operation. Real HAO supports
are smaller than full AO supports, so even the projector speedup differs.

Aggregate measured seconds:

| Phase | 241/block-LBFGS old → new | FeCl2/block-LBFGS old → new | FeCl2/NEO old → new |
|---|---:|---:|---:|
| Projection, including shared setup | 0.016255 → 0.009020 | 2.781550 → 0.512382 | 0.053793 → 0.009551 |
| Normalized curvature | 0.001158 → 0.000928 | 1.415959 → 1.408163 | 0.027099 → 0.027181 |
| Positive block factorization | 0.010176 → 0.008212 | 13.869701 → 13.805660 | 0.267530 → 0.266309 |
| Whole chart, inclusive | 0.062277 → 0.045973 | 38.292064 → 35.822469 | 1.047656 → 1.000728 |
| Trial forward | 0.483378 → 0.476068 | 18.570598 → 18.503765 | 0.117162 → 0.128618 |
| Trial gradient | 0.119953 → 0.115025 | 3.668814 → 1.426938 | 0.022300 → 0.022064 |

Chart totals include their component phases and must not be added to them.
Likewise block eigensolution is included in block factorization; for FeCl2/
block-LBFGS it alone is 13.001516 → 12.962213 s. These block factorizations were
not optimized. Another approximately 20 s of chart work is outside the above
three components; the current timers do not subdivide it, so do not attribute
that residual to a particular chart operation without further measurement.
For NEO, the majority of total time is outside the chart and trial phases;
this audit does not break down that remainder into HVP/Krylov/response costs.

For FeCl2/NEO, the original projection was just about 0.46% of total time.
Saving 44 ms cannot turn an 11.6-s calculation into a much faster one. For
FeCl2/block-LBFGS, projection saves 2.269 s and trial gradients save 2.242 s,
together explaining most of the measured 4.89-s wall reduction.

FeCl2/block-LBFGS Armijo counters: accepted 256 → 255, rejected 402 → 403.
Total trial gradient calls drop from 659 to 256. There is one additional
completed trial outside the Armijo accepted-count scope. For 241/block-LBFGS,
there are no rejected Armijo trials; energy-only screening cannot save any
rejected-trial gradient work there. NEO does not use this Armijo path.

## Memory: live heap versus resident pages

Independent heaptrack captures, kept separate from native timing/RSS, yielded
exact peak live heap bytes via Massif-compatible reports:

| Capture | Old peak bytes | New peak bytes | MiB |
|---|---:|---:|---:|
| 241/block-LBFGS | 1,043,838,599 | 1,043,838,599 | 995.482062 |
| FeCl2/NEO | 425,359,291 | 425,359,291 | 405.654231 |

Peak allocation stacks contain no OrbitalChart allocations at the global
peak; large integral/input-loading allocations precede the optimizer. Both
versions report the same 31.14 kB unfreed total at process exit, not a new
projector leak. Do not use merged allocation-site peak totals: heaptrack's
help explicitly warns that merged peak consumption is not correct. Exact
global peaks above come from the Massif byte series, not summed site peaks.

Largest owned projector matrix payload is only 19,872 bytes for 241 and
328,136 bytes for FeCl2. Referenced input matrices are not counted as owned.
These payloads cannot explain the previous approximately 10–17 MiB RSS increase.

Allocator ablation (`MALLOC_ARENA_MAX=1`, separate paired nodes):

| Case | Default RSS MiB, old → new | One-arena RSS MiB, old → new | One-arena wall s, old → new |
|---|---:|---:|---:|
| 241/block-LBFGS | 962.1 → 969.6 | 772.3 → 773.3 | 3.25 → 3.35 |
| FeCl2/NEO | 449.2 → 448.4 | 350.4 → 355.5 | 10.25 → 10.18 |

This demonstrates substantial allocator-sensitive resident memory, consistent
with retained/free arena pages and parallel allocation variability rather
than added persistent numerical data. It does not attribute every residual
MiB difference to a particular arena. 241 slows substantially with one arena,
so **do not adopt this environment setting as an application default**.

## Conclusion

Projection reuse and rejected-gradient screening reduce the intended work.
They are small parts of these total workloads; the isolated benchmark was
not a basis for expecting a sevenfold overall improvement. Same-node native
controls do not reproduce the earlier slowdown, and exact peak live heap is
unchanged. Substantial RSS variability is allocator-sensitive. The remaining
large chart costs, block decompositions, trial forwards, and NEO model-solving
costs require separate attribution before selecting another optimization.
