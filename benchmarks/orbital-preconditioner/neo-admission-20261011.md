# Historical coupled-model admission ablation

## Controlled change

Starting from production revision `0623d1c`, the experimental patch changes
only the model selection in the NEO backend: use coupled NEO whenever the
structure tangent space is nonempty, as in `c2776c4`. It does not restore
historical geometry, derivatives, preconditioning, response accuracy or final
convergence checks. The production source and executable remain unchanged.

`neo-admission.sbatch` builds the patched executable on a compute node and
restores the execution mirror before either run. Both native variants run
sequentially in the same allocation, with the same compact affinity, 32 OpenMP
threads, one BLAS thread, exact Libcint integrals, Davidson, and the historical
NEGPIA input. Projected-gradient and adjacent-energy thresholds remain
1e-3 and 1e-7 Eh. There is no added optimizer iteration limit.

Slurm job `281618` completed on Hanhai25 `anode042`; both application runs
exit zero and explicitly report convergence. This is one matched pair, not
a repeated timing study. GNU time supplies native process peak RSS; Slurm
CPU accounting is incomplete and is not used to infer CPU utilization.

## Results

| Measurement | Current admission (core-only) | Historical admission (coupled) |
| --- | ---: | ---: |
| Accepted points | 17 | 6 |
| Orbital HVPs | 168 | 45 |
| Structure actions, including iterative certification | 0 | 293 |
| Rejected trials | 4 | 0 |
| Total wall time, s | 19.452345 | 8.426375 |
| SCF wall time, s | 17.626676 | 6.646346 |
| Peak RSS, MiB | 3765.51 | 3726.29 |
| Initial energy, Eh | -435.698724650239 | -435.698724650239 |
| Final energy, Eh | -436.262505347363 | -436.262507052760 |
| Final projected gradient infinity norm | 2.71668521e-4 | 5.38250561e-4 |
| Final accepted energy decrease, Eh | 6.41070983e-8 | 7.10393806e-8 |

The coupled variant is 2.31 times faster end-to-end and uses 73.2% fewer
orbital HVPs in this pair, despite adding structure actions. It recovers the
historical six-point convergence without reverting correctness fixes.
Both final projected gradients and energy changes satisfy the unchanged
outer criteria. Their final energies differ by 1.705397e-6 Eh; equal stopping
criteria do not imply identical final energies, and no tighter endpoint
agreement is claimed.

The coupled first-step actual decrease is 0.474794220168 Eh, close to the
historical 0.474794220727 Eh. The production core-only decrease is
0.474511841843 Eh. This strengthens the attribution to the changed model,
not merely to the final stopping criterion.

## Conclusion and limits

For NEGPIA, switching from coupled to core-only to reduce estimated action
cost caused a substantial convergence regression. This switch changes both
the local-active orbital Hessian contribution and the orbital/structure
coupling: `apply_orbital_coupling` enables `local_active_response`, whereas
the core-only branch disables it. This experiment therefore does not isolate
structure coupling as the cause. Model admission must account
for model quality and total nonlinear cost, not only per-action dimensions.
This does not establish that complete coupling is economical for every
large active space. Production admission has not been changed: other
core-only systems must be tested before choosing a general replacement.

The next experiment should compare the same two models on the remaining
core-only panel systems, retaining all current tolerances and memory controls.
Avoid adding a NEGPIA-specific exception or reinstating old loose residual
certification.

Full hashes and numerical summaries are in
`results/history/neo_admission_20261011.csv`. Raw logs are saved locally under
`results/studies/locality-20261011/neo-admission-281618/` and remotely under
`/home/guqqgroup/taoxia/xmvb-runs/native-20261011/neo-admission-281618/`.
