# Sigma and directional-HVP audit

## Outcome

This round did **not** reproduce the large exterior-transform speedup of the
previous round. Five converged native runs show mostly small timing changes,
including a QOCGUM regression. Peak RSS does not consistently decrease.
Do not interpret these single-run differences as an established speedup.

Retained changes:

1. Sigma's beta reduction visits contiguous rows within each worker tile,
   instead of crossing all worker buffers for each scalar output. The first
   existing worker tile holds the sum. Per-element worker addition order is
   unchanged; output ownership and the barriers between tiles are unchanged.
   No new buffer or thread-private full spin-string carrier is allocated.
2. Directional H/S projection applies the beta transpose transform directly to
   non-owning column blocks. Two carrier-block allocations, copies and copybacks
   are removed from this path. Transform order and formulas are unchanged.
3. Separate optional timers identify directional exterior transforms, integral
   adjoints and one-body-generator adjoints.

No optimizer policy, threshold, HVP formula, or outer convergence condition was
changed. The existing sigma workspace budget was not increased. These edits do
not establish a global memory bound for every other code path.

## Rejected experiment

Build 1250768 grouped beta strings into larger GEMMs within the existing
workspace budget. Array 1250771 converged, but gave 2.301, 8.122, 14.660,
20.010 and 34.944 seconds for 241, FeCl2, QOCGUM, CERRAS and LOFLEA.
It showed no stable benefit and was removed, including its extra helper.
There is no retained switch or alternative batching path.

## Validation and native results

Build 1250784 passed eight test launches in Slurm: direct-CI at 1 and 32 threads,
Davidson, keyframe, CIAH, coupled structure, sparse quotient and Hessian diagonal.
The tile-boundary test now uses 252 spin strings and two RHSs, covering a partial
alpha tile, uneven beta ownership, dynamic and nested OpenMP teams. Existing
tests cover integral-adjoint bilinear identities and directional derivatives.

Native array 1250785 used Hanhai22, 32 bound OpenMP cores, single-thread BLAS,
exact libcint integrals, Davidson and analytic-diagonal NEO. Input and basis
SHA-256 hashes match saved array 1250728, the previous retained version.
The old executable was not rerun. Each case landed on the same node name as
its saved reference; shared-node conditions still prevent attributing small
timing differences entirely to this patch.

Wall time is the application's end-to-end time. RSS is GNU time process peak,
not Slurm requested memory. Values are before / after.

| Case | Steps | HVP | Wall / s | RSS / MiB | Energy change / Eh |
|---|---:|---:|---:|---:|---:|
| 241 | 3 / 3 | 21 / 21 | 2.122 / 2.069 | 882.3 / 820.7 | 0 |
| FeCl2 | 12 / 12 | 251 / 256 | 8.133 / 7.963 | 491.9 / 491.2 | -6.23e-10 |
| QOCGUM | 4 / 4 | 52 / 52 | 12.632 / 13.165 | 2030.8 / 1915.4 | +1.02e-12 |
| CERRAS | 5 / 5 | 51 / 51 | 18.341 / 18.174 | 1928.3 / 1957.1 | -9.66e-13 |
| LOFLEA | 4 / 4 | 34 / 34 | 34.875 / 34.131 | 2745.7 / 2874.6 | -1.02e-12 |

FeCl2 uses five more HVPs in this run; this patch does not guarantee identical
optimization trajectories. Final energies agree well within the unchanged
energy tolerance. There is no demonstrated lower-RSS result for LOFLEA.

`sigma-comparison.csv` includes final energies, node names and input checks,
including the rejected experiment. Regenerate its contents with:

```sh
python3 compare-sigma.py sigma-batched-1250771 sigma-reduction-1250785
```

## Actual directional call chain

The two-macro LOFLEA diagnostic 1250786 contains 13 HVP expansions:

| Scope | Calls | Inclusive wall / s |
|---|---:|---:|
| Trial forward | 3 | 9.006 |
| HVP apply | 13 | 2.959 |
| Sigma | 366 | 9.006 |
| Directional exterior transform | 0 | 0 |
| Integral adjoint | 32 | 2.701 |
| Generator adjoint | 58 | 0.358 |

These scopes overlap and must not be summed. The short diagnostic stops at its
explicit two-macro limit; it is not a failed full convergence test. Profiling
and native LOFLEA shared an allocated node, not CPU cores. No sampled VTune
profile was obtained under the cluster's current ptrace policy.

The main HVP obtains orthogonal coefficient and sigma directions through
one-body generators and reuses them in backward. The directional shear API is
not invoked in this measured trajectory. Rewriting that API would miss the
observed hotspot. The optional path still has callers, so it is not dead code.

The next worthwhile design question is whether the response sigma and integral
adjoint can share their one-body scatter within a bounded tile. Both consume
the same response coefficient carrier in `active_integral_response_adjoint`.
That requires a fused action/adjoint interface and a memory analysis; no savings
from it are claimed or implemented in this patch.
