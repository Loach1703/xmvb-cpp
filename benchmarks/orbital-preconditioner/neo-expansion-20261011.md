# MnF2 NEO expansion audit

## Scope

Base revision: `808451d`, plus scalar diagnostic instrumentation only.
Hanhai25 job `281305`, anode020, one 32-core allocation, compact affinity,
one BLAS thread, exact Libcint, Davidson, one-electron preconditioner.
Gradient and energy thresholds remain 1e-3 and 1e-7 Eh. Existing Release
build flags are unchanged. No extra coupling, metric, or response actions
are introduced by `XMVB_NEO_TRACE=1`; all recorded quantities already exist.
Thread-local diagnostic context holds two integers and a source label only.

The same binary runs once without scalar logging and once with it, in that
order. Both use the existing macro `--neo-trace` callback. These are single
runs, not statistically established timing differences or identical
intermediate trajectories. Tests `test_neo`, `test_neo_response_solver`, and
`test_shifted_metric_preconditioner` passed on the compute node.

| Scalar diagnostics | Steps | Orbital HVPs | Final energy, Eh | Application wall, s | Process peak RSS, MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| Off | 14 | 669 | -1348.893353119315 | 12.681985 | 515.1 |
| On | 14 | 669 | -1348.893353119319 | 12.729823 | 499.1 |

RSS differences are not a memory optimization: instrumentation adds no
large physical arrays, and resident high-water marks vary between runs.

## Findings

1. All 669 orbital basis proposals are independent: 14 initial directions
   and 655 preconditioned orbital-residual directions. There are no rejected
   dependent proposals or canonical-coordinate extensions. The minimum
   retained Euclidean norm fraction is approximately 0.089. This excludes
   rejected/near-zero proposals as the explanation for this run, but does
   not measure the physical conditioning of the retained subspace.
2. There are 25 candidate-local response checks across 18 radius trials.
   Seven checks trigger structure-space enrichment. Hence this mechanism
   affects the trajectory but is not hundreds of repeated enrichments.
   These check counts are not total iterative structure-operator actions.
3. Steps 10 through 12 use 137, 231, and 102 orbital HVPs: 470 of 669,
   about 70.3%. Steps 10 and 11 are boundary steps. Their shifted
   preconditioner is invoked 69 and 127 times respectively, with zero inner
   PCG iterations in every invocation. The initial shifted block inverse
   already satisfies the current loose direction-generation budget. Step
   12 is unshifted and directly uses the block inverse. This confirms that
   the expensive part is many outer orbital extensions, not a nested
   accurately solved preconditioner in these three steps.
4. The first radius-1 trial needs 51 HVPs before discovering a predicted
   decrease of 27.6261 Eh corresponds to an actual increase of 0.7780 Eh.
   Radius 0.25 is also rejected; radius 0.0625 is accepted. The latter two
   reuse the 51-direction space and require no additional orbital HVPs.
   The retained model is reused correctly; its initial nonlinear validity
   is poor at the proposed radius.
5. Step 7 rejects radius 0.5 after 32 HVPs and radius 0.125 after another
   six; radius 0.03125 is accepted after six more. Four rejected trials
   collectively precede 89 newly generated HVPs. Do not call all 89 wasted:
   their basis/images remain available for the accepted smaller-radius
   candidates. Per-attempt logging separates this from restarting a model.

## Interpretation and next measurement

The evidence supports inspecting the preconditioned relaxed curvature and
the cost of obtaining a boundary candidate before nonlinear validation.
It does not support changing SVD drivers, canonical-direction handling, or
eliminating candidate-local response certification as the main remedy.

A one-electron block inverse does not explicitly represent two-electron
inter-orbital curvature or the structure Schur correction. The measured
direction-generation budget accepts that inverse without refinement in
the expensive late steps. This is a concrete limitation to investigate,
not proof of a particular condition number or that merely tightening PCG
will improve total time. Even an exact inverse of the local model would
still omit its missing relaxed curvature.

Preserve the physical KKT and final stationarity certificates. Before
changing the solver, compare the projected physical residual and model
decrease obtained per added direction, using small projected spectral
information and existing images. A cheap nonlinear boundary trial needs
a mathematically stated sufficient-decrease contract, not a new arbitrary
HVP cap or weaker final convergence threshold.

## Reproduction

Submit `neo-expansion.sbatch SOURCE_DIR` on Hanhai25. Parse scalar records:

```sh
python3 benchmarks/orbital-preconditioner/neo_trace_results.py \
  results/studies/locality-20261011/neo-expansion-281305/trace/run.err \
  --csv results/studies/locality-20261011/neo-expansion-281305/trials.csv
```

Local evidence includes native output, resources, affinity, input/binary
hashes, macro traces, and compute-node test logs under that run directory.
Remote evidence is under
`/home/guqqgroup/taoxia/xmvb-runs/native-20261011/neo-expansion-281305/`.
