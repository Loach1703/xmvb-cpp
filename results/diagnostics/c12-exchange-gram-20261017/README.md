# Exact exchange Gram correction: numerically qualified, same-node timing pending

## Change and mathematical necessity

The response integral adjoint needs A^T D B, where D weights each off-diagonal
spin-exchange representative by two and each diagonal representative by one.
Use the exact identity **A^T D B = 2 A^T B - A_diag^T B_diag**. The derivation,
arbitrary-selected-space scope, cost and numerical limitations are recorded
in `article/exchange_gram.md`.

Diagonal indices are gathered during existing target compaction. A shared
O(number_of_groups) integer list replaces no physical amplitude and does not
depend on a residual or coefficient threshold. Rank-one diagonal corrections
are applied to the existing small partial Grams or disjoint output-channel
blocks. There is no new numerical slab or per-worker matrix allocation.

The common accepted-state channel slab can therefore remain unweighted across
all temporal requests, until its final in-place directional image consumes it.
Full-slab weighting and later unweighting (or unweighting the directional image)
are removed. Output exchange reconstruction, physical orbital and CI forcing,
integral adjoints, precision, workspace budgets, tile planning, OpenMP
ownership, convergence criteria and optimizer policy remain unchanged.

This is an exact alternative contraction, not an asymptotic reduction of the
complete sigma/HVP. Removing O(N P) passes introduces O(d P^2) diagonal work;
it need not be faster for small, diagonal-heavy tiles. Performance qualification
is separate from proving the identity or passing numerical tests.

## Workload and numerical qualification

Hanhai22 job **1264878**, 32 physical cores, GCC14.1/C++23 Release production
flags, BLAS threads=1. The compiled `sigma_tiles.cpp` hash matches the local
candidate. No unrelated local RI changes were transferred.

The C12 workflow uses CAS(12,12)/COV132 operator fixtures plus complete
**LOFLEA HAO** convergence: exact Libcint/cc-pVDZ, 150 basis functions, 132
structures and 924 strings per spin, ISCF=9/Davidson. This is not a separate
carbon-C12 molecular input. Original orbitals, thresholds and ITMAX=2000.

Input SHA256:
`092d545498678b0cf772c5c84b9eab674b694375558cd391a1b739b7bdcc35eb`.
Basis SHA256:
`d53e7e02d55fb99cdf823d50cf13fec375e0cbe1220ff61816aee7ab38e38ec9`.

Job 1264878 completed exit0; **9/9 regressions passed in 3.48 s**, including
both Gram partition plans, different exchange parities, unequal spin sectors,
multiple RHSs, temporal directions and shuffled/compact sparse output.
Independent CAS12 checks retain:

- HVP/gradient finite difference error: 2.5933568576767072e-11.
- Lagrangian finite difference error: 6.9466100081384818e-12.
- Physical structure-response image error: 3.4002682088792616e-16.
- Structure-response integral-adjoint error: 1.2707391900175511e-15.

The optional complete-carrier jet replay is not enabled; its placeholder zeros
are not additional evidence. Finite differences and structure-response
reference checks are executed independently at their existing tolerances.

## Native results: different nodes, no attributed speedup

| Version | Job | Node / affinity | Steps / HVP / responses | Final energy / Eh | SCF / s | Application / s | Peak RSS / MiB |
|---|---:|---|---:|---:|---:|---:|---:|
| Saved parent 505e5405 | 1264871 | anode40 / 60-91 | 4 / 36 / 36 | -399.106646267329 | 19.845801 | 22.626016 | 2687.16 |
| Exchange Gram candidate | 1264878 | anode42 / 96-127 | 4 / 36 / 36 | -399.106646267331 | 16.871700 | 19.492537 | 2719.13 |

Both nodes report AMD EPYC 9654 / two 96-core sockets / eight NUMA nodes.
Node load and CPU/NUMA allocation differ. Several **unchanged** phases also
accelerate. These timings do not establish a same-node algorithmic speedup or
RSS reduction. No old binary was executed to obtain a new baseline.

Separate phase runs preserve 129 sweeps, 1171 batches, 1120445436 logical source
elements and 1068 product projections. Coupled-action scope changes from
6.522412 to 5.636303 s, Gram from 0.510799 to 0.410206 s, and unchanged AO-H1e
scope from 2.021217 to 1.799605 s. These are nested scopes,
not independent additive wall-time shares. The cross-node improvement must
not be attributed entirely to the exchange change.

## Scheduling and frozen same-node follow-up

The original requested node anode40 became fully allocated. The same pending
job 1264878 was updated, without restarting it, first to anode10 and then
anode42. Anode10's apparently idle cores belonged to an exclusive owner; no
workload ran there. Anode42 was scheduler-accessible and ran the job normally.
No time limit was shortened and no other user's jobs were changed.

**Job 1264881** is queued on **anode40**, using `c12_regression.sbatch --bundle`
with the exact frozen 1264878 binary pair, source manifest and original input
files. This mode skips rebuilding mutable future source; native and phase
calculations plus CAS12 checks use the same qualified candidate. It is a new
version run, not an old-baseline rerun. The existing default build-and-test and
tests-only modes are preserved. Bash syntax validation passed locally; actual
bundle-mode execution remains pending with the Slurm job.

Raw logs, build cache, hashes, affinity and numerical outputs are archived,
excluding binaries. `comparison.csv` uses the existing native parser and
`phase-summary.json` the existing phase parser. Check Slurm job 1264881 rather
than assuming it completed or resubmitting because an observation timed out.

Follow-up: job 1264881 was cancelled while still pending, before executing any
calculation, because the workspace-lease candidate superseded it. Frozen
latest-candidate job **1264892** subsequently completed on anode40; see
`../c12-workspace-lease-20261017/README.md`. The anode42 data above remain
unchanged and no old baseline was rerun.
