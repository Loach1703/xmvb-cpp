# Unclipped CIAH tail: controlled ablations

## Decision

Retain the optimization policy from `146dfc87`. None of the three tested
changes gives a satisfactory improvement across the five cases. Both source
experiments were removed, not retained behind options. The existing
one-electron preconditioner remains an existing alternative, not the default.

These are converged native runs on Hanhai25, 32 physical cores, BLAS threads 1,
exact Libcint integrals, CIAH/Davidson. The input files, outer tolerances,
coupled physical HVP, nonorthogonal coordinates and CI transport are unchanged.
The one-electron test changes only the existing preconditioner selection;
the other two use the analytic Hessian diagonal. No old run was repeated.
Builds and calculations ran in Slurm compute jobs, not locally or on login nodes.

| Variant | 240 HVP | 241 HVP | FeCl2 HVP | F2 HVP | MnF2 HVP |
| --- | ---: | ---: | ---: | ---: | ---: |
| Retained, job 291714 | 27 | 29 | 170 | 18 | 160 |
| One-electron, job 291729 | 24 | 38 | 145* | 16 | 135 |
| Current-gradient seed, job 291751 | 25 | 29 | 162 | 16 | 218 |
| Base inverse only, job 291758 | 27 | 31 | 172 | 23 | 182 |

*FeCl2 stops at `-2181.617645246294 Eh`, **6.819625 mEh higher** than the
retained solution `-2181.624464871226 Eh`. Its smaller action count is not a
successful replacement. All other endpoints differ from their retained
counterparts by less than `1.3e-9 Eh`.

| Variant | FeCl2 steps / wall / RSS | MnF2 steps / wall / RSS |
| --- | --- | --- |
| Retained | 13 / 3.68 s / 444.67 MiB | 11 / 2.66 s / 506.48 MiB |
| One-electron | 6* / 2.76 s / 431.08 MiB | 9 / 2.43 s / 521.63 MiB |
| Current-gradient seed | 12 / 3.40 s / 433.69 MiB | 13 / 3.32 s / 513.48 MiB |
| Base inverse only | 14 / 3.77 s / 434.86 MiB | 11 / 2.75 s / 506.81 MiB |

These are single observations on different nodes with the same CPU family,
not repeatability-qualified timing ratios. The action counts and endpoints
already reject the proposed defaults; small wall differences are not
credited as speedups. Raw paths, hashes, node names, full energies and all
five cases are in `summary.csv`.

## What was changed

- Current-gradient seed: replace the transported previous macro seed by the
  current preconditioned gradient. Original coarse correction remains active.
- Base inverse only: restore the original seed and bypass only the auxiliary
  coarse balancing inverse. Full orbital/CI coupling remains in every HVP.
  Sample observation still runs in this diagnostic, so its timing does not
  measure the savings of removing sample storage/observation altogether.

The current-gradient-seed and base-inverse-only sources passed the seven
focused tests in build jobs 291750 and 291757, respectively. The compiled
experimental executables are identified by hashes in the CSV; they are not
the retained production policy.

## Tail diagnosis

The retained FeCl2 last three macros use 25, 27, 29 HVPs; MnF2's last four
use 18, 20, 22, 26. These are unclipped microsteps, so replay cannot remove
their HVP cost. Late FeCl2 logs show projected Galerkin residuals near
`1e-18` while the full residual is of order `1e-4`. This is a missing
complement direction, not an inaccurate small projected solve. See the
derivation in `article/ciah_clipped_path.md`.

The old diagnostic named `coordinate_cosine` divided a Euclidean dot product
by a coordinate norm and a Euclidean norm. It is not a cosine and must not
be used to diagnose angular quality. The diagnostic is corrected to a true
`euclidean_cosine`; this affects trace output only, not the optimization.

After restoring the retained policy and making this trace-only correction,
build job **291788** passed all seven focused tests on Hanhai25 (exit 0).
Its saved build/test artifacts are under
`results/diagnostics/ciah-tail-20261017/retained-build/`. The retained native
baseline was not rerun.

No large-system improvement follows from these tests. The remaining target
is more effective complement directions; neither removal of coarse
corrections nor replacement of the macro seed has established that benefit.

## Reproduce the table without computing

Run `python3 benchmarks/ciah-tail-20261017/summarize.py`. It reads the saved
local native logs, checks exit/convergence status and matching input hashes,
and prints CSV. `ritz_solves` counts AH solves, not physical microsteps.
The runner in `../ciah-path-20261017/run.sbatch` accepts an optional existing
preconditioner name for controlled experiments; no new production switch
was added.
