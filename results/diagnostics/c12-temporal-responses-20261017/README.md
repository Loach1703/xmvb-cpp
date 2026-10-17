# Temporal orbital/CI response sharing

## Scope and invariant

Production change relative to `fe59d7a6`: a CIAH orbital direction and its CI
direction share a generated graph, accepted read/scatter and caller-owned
exterior-transform panel caches. They retain separate Hamiltonian images,
forcing, metric images and integral/generator cotangents. State weights,
normalization/gauge terms and orbital adjoint pullbacks are unchanged.
No convergence or residual threshold, model policy, coordinate formula or
physical structure response was removed. No global spin-product table or
thread-private carrier was added. The 64 MiB sigma workspace bound is unchanged;
the existing exterior panel caches have their own unchanged budgets.

See [theory](../../../article/temporal_ci_responses.md) for the lifetime proof.
The in-place two-slab schedule applies to the ordered CI-then-orbital sequence,
not to standalone directional actions. The latter retain ordinary GEMM.

## Native end-to-end results

LOFLEA, exact Libcint/cc-pVDZ, HAO, CAS(12,12), COV 132 structures,
924 alpha and 924 beta strings, CIAH/Davidson, original guess. Hanhai22 anode40,
32 physical cores, CPU affinity 60-91, GCC14.1.0, C++23, Release
`-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`, Eigen BLAS/LAPACKE enabled,
OpenMP32 and BLAS1. Original ITMAX2000 and convergence thresholds retained.
Times include complete convergence, not a single action or SCF step.

| Version / saved job | Steps | HVP | Final energy / Eh | SCF / s | Application / s | Peak RSS / MiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Scatter baseline / 1264636 | 4 | 36 | -399.106646267331 | 28.207024 | 30.863387 | 2824.68 |
| Halo baseline / 1264676 | 4 | 36 | -399.106646267332 | 28.794006 | 31.677382 | 2730.44 |
| Saved repeat of same halo binary / 1264685 | 4 | 36 | -399.106646267330 | 27.658785 | 30.361555 | 2775.21 |
| Shared physical chain / 1264713 | 4 | 36 | -399.106646267330 | 26.043647 | 29.009265 | 2725.68 |

The native application is 6.0% faster than the saved scatter baseline and
4.5-8.4% faster than the saved halo measurements. This is a single new native
measurement, not a statistically qualified speedup interval. RSS comes from
GNU time (KiB divided by1024), not Slurm's sampled MaxRSS. No substantial RSS
reduction is claimed. Input and basis hashes match the saved runs.
Every native run retained four accepted steps, 36 HVPs and 36 physical CI
response actions. Their per-step HVP counts remain 6,17,10,3.

## Phase-instrumented attribution

These are separate instrumented convergence runs, not the native timings above.
Scopes overlap and cannot be summed. The orbital HVP scope now includes the
shared CI sweep: its old and new boundaries are different. Compare the full
`ciah_images` scope rather than treating the moved work as an orbital slowdown.

| Quantity | Saved 1264676 | New 1264713 |
| --- | ---: | ---: |
| Application wall / s | 30.913943 | 28.592445 |
| Coupled image wall / s | 11.990040460 | 9.768214716 |
| CI response kernel invocations | 88 | 52 |
| Generated sweeps | 165 | 129 |
| Generated batches | 1639 | 1171 |
| Source elements read | 1420939020 | 1120445436 |
| Tile preparation wall / s | 1.331151489 | 1.035013168 |
| Tile read wall / s | 2.644525129 | 1.464799599 |
| Tile compute wall / s | 6.374195725 | 5.584977433 |
| Tile write wall / s | 2.129010091 | 2.114584340 |

The physical response count did not decrease: 36 duplicate kernel invocations
were merged. The full coupled-image scope improved18.5%, source traffic21.1%
and generated batches28.6%. These are traced contraction counters and scoped
wall times, not profiler sampling percentages. Existing profiler availability
limitations are unchanged; no new CPU sample or allocation capture is claimed.

## Operator tests and intermediate designs

`run-1264713/ctest.out`:9/9 passed. The direct-CI tests cover parity0/+/-1,
unequal spins, multiple RHSs, strict support, forcing and adjoint-only consumers,
tiny halos and nested teams. Product-jet tests compare mixed orbital/CI
responses for cached and uncached accepted images, nonclosed selected spaces,
multiple states, zero weights and exact-zero directions. The CAS12 numerical
operator differences are below1e-15.

Job1264750 adds a complete physical-chain comparison of mixed, pure-CI and
pure-orbital directions, including tangent Hessian/metric, gauge multipliers,
orbital forcing and final adjoint pullback. It passed at the unchanged3e-11
test tolerance. This is equivalence qualification against existing independent
separate calls, not independent absolute derivative validation.

`--temporal12` is a warm current-API microbenchmark with dense **test-only**
924x924 oracle arrays. It is not the low-memory production caller. Job1264695's
three-slab sequence took0.3072s versus0.2777s separately. The ordered two-slab
sequence in1264696 took0.2061s versus0.2568s, but in1264713 took0.4126s versus
0.2921s. The microbenchmark gain is therefore not robust and is not used as an
end-to-end speedup claim. Deterministically, accepted read callbacks fell22→9;
all image/cotangent differences stayed below1e-15. Native complete-convergence
timing and the phase counters above qualify the integrated change.

The physical chain was still separate in1264695/1264696; their full-run results
are included in `comparison.csv`, not presented as production fusion results.
Job1264700 failed to build due to one internal call missing the new explicit
null optional request; it executed no numerical workload. The call was fixed
before1264713. A failed `sbatch --wrap` submission for the extra test created no
job;1264750 used the archived script instead.

## Reproduction and provenance

Only current candidates were executed. All old measurements above are saved
artifacts; no old executable was rebuilt or rerun. Three unrelated local dirty
files (`.gitignore`, RI AO operator and its unit test) were not transferred or
included in this change. The compiled changed production/test files all match
1264713's source hashes; the additional coupled test has1264750's own hash.

The application binary SHA256 is
`21574a3a0dfbf6baabbe9246836f04bd034f7bbb2fa095accf970d61461f5f11`.
The input SHA256 is
`092d545498678b0cf772c5c84b9eab674b694375558cd391a1b739b7bdcc35eb`;
the basis SHA256 is
`d53e7e02d55fb99cdf823d50cf13fec375e0cbe1220ff61816aee7ab38e38ec9`.
The run artifacts retain CMake flags, source/binary hashes, node/CPU/module/
affinity environment, input, output, convergence trace, native RSS and test logs.
Frozen binaries remain on Hanhai22 under
`/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015/c12-regression-1264713/bin/`;
they are intentionally excluded from the repository archive.

Unpack `evidence.tar.gz` here to inspect the raw new runs. `summarize.py` produces
`comparison.csv` from those and the existing scatter/halo artifact directories.
`tests/mcscf/summarize_c12_cost.py` produces `phase-summary.json` without running
any calculations. `tests/mcscf/c12_regression.sbatch` is the convergence launcher;
`test-coupling.sbatch` launches only the extra full-chain unit test.

## Remaining necessary work

Each independent response still needs its Gram product, Hamiltonian image and
selected-space adjoint projection. Generator-adjoint consumers still traverse
their direction-specific arguments separately; that is the next mathematical
fusion boundary to examine. No claim of an asymptotic-order reduction or a
finished global performance optimization is made by this change.
