# Certified occupation-block metric inverse in CIAH

## Mathematical change

The factorized complete-space structure overlap inverse already existed in
`structures/orthogonal_ci/metric.cpp`. CIAH's horizontal dual norm nevertheless
preconditioned its Riesz CG with the raw diagonal. Full variational spin-sector
completeness, not carrier completeness, is needed to reuse the factorization.

With independent tangent basis $Z$, selected roots $C$ and overlap $M$,

$$
Z^TMC=0,\qquad G_s=Z^TMZ,\qquad
G_s^{-1}=Z^T\left[M^{-1}-C(C^TMC)^{-1}C^T\right]Z.
$$

The subtraction is essential: simply projecting $M^{-1}$ is not the exact
restricted inverse. The derivation and complete-space exterior factorization
are in `article/neo_structure_metric.md`. Selected clusters use all roots in
this correction, with no additional weight factor because the coordinates
already include $\sqrt{2w_s}$.

The accepted point now owns the admitted occupation-pattern factors once.
CIAH and applicable response consumers borrow this same immutable object.
The factor/workspace cap remains the existing 64 MiB contract. There is no
global structure overlap or Hessian allocation and no per-thread copy.
Noncomplete or nonadmitted spaces continue to solve their own complete metric
matrix-free; their diagonal model is only a preconditioner.

This replaces a preconditioning model, not the physical metric or its norm.
The same solver still explicitly certifies the true residual at
sqrt(machine epsilon). No tolerance, microstep budget, AH equation, Hessian
action, coordinate chart, retraction or outer stopping rule changed.

## FeCl2 measurement

Hanhai22 32 bound cores, original guess, exact def2-TZVP, Davidson,
analytic diagonal, gradient 1e-3 / energy 1e-7 Eh. Both native runs were on
anode31. The older run is reused; no old executable was rerun.

| Quantity | Archived 1249395 | New 1249408 |
|---|---:|---:|
| Macros | 14 | 14 |
| Orbital HVP | 281 | 280 |
| Energy / Eh | -2181.624464869471 | -2181.624464869629 |
| GNU wall / s | 10.55 | 9.81 |
| GNU RSS / MiB | 512.3 | 496.6 |
| Complete structure norm calls | 505 | 504 |
| Norm wall / s | 0.82504 | 0.11507 |
| Overlap actions in norm solves | 22742 | 924 |
| Largest CG iteration count | 53 | 1 |
| Largest certified relative residual | 1.4882e-8 | 2.3700e-15 |

The norm scope is 7.17 times faster and uses 24.61 times fewer overlap actions.
Total wall time decreased 7.0% in these single same-node observations; no
statistical timing claim is made. The difference of one HVP is not a material
convergence improvement. Peak RSS did not increase, but its observed decrease
must not be attributed to the 7.7 KiB factors without allocator evidence.

FeCl2 factors total 7,840 bytes, 308 spatial occupation blocks, largest block
5. The endpoint differs from the retained PySCF coupled-Newton reference
-2181.624464883733 Eh by +1.4104e-8 Eh. The microstep-budget tail remains;
this optimization does not complete the convergence-efficiency goal.

All five compute-node unit executables passed. New tests independently
assemble small one-state and equal-weight two-state coordinate metrics,
verify the exact inverse and norm, check sharing and subquadratic factor
storage, and retain true-residual certification. A context without variational
completeness cannot admit the factors, even when its carrier is complete.

## Small-system regression

Frozen 1249408 executable, parallel Slurm array 1249412. All converge with
exit zero; archived run 1249320 supplies the comparison without rerunning it.

| Case | Archived macro/HVP | New macro/HVP | New energy / Eh | New GNU wall / s | New GNU RSS / MiB |
|---|---:|---:|---:|---:|---:|
| F2 | 3 / 23 | 3 / 23 | -198.751155830526 | 0.33 | 17.5 |
| 241 | 3 / 21 | 3 / 21 | -230.720590392836 | 2.76 | 829.1 |
| MnF2 | 10 / 144 | 10 / 144 | -1348.893353223967 | 5.00 | 496.6 |

Energy differences are at approximately 1e-12 Eh or below. All nonzero
structure Riesz solves take one iteration. Factor storage is 24, 2,040 and
3,648 bytes respectively. These cases ran on different nodes from archived
measurements; do not infer a universal timing improvement from this table.

## Bounded iron-porphine follow-up

Array 1249418 uses the frozen new executable and the original input orbitals:
singlet and triplet separately, four macros each. This is not an optimized
restart, not a convergence result, and not a new long run. Its purpose is to
measure factor admission, certified norm costs and RSS at 32,670 / 50,820
structure dimensions. Final energies and full convergence efficiency still
require subsequent demonstrated direction/model improvements.

Both tasks have now terminated at the intended four-macro diagnostic cap
(exit 2, not converged):

| Spin | Macros / HVP | Energy / Eh | GNU wall / s | GNU RSS / GiB | Structure norm / s | HVP expansion scope / s |
|---|---:|---:|---:|---:|---:|---:|
| Singlet | 4 / 36 | -2244.717977952240 | 424.92 | 38.621 | 2.139 | 19.257 |
| Triplet | 4 / 42 | -2244.769462512970 | 621.89 | 38.630 | 2.100 | 26.289 |

All nonzero structure metric solves need at most one CG iteration, with
maximum certified relative residual 2.14e-15 / 2.67e-14 respectively.
Factor storage is 982,080 / 3,067,680 bytes, shared at the accepted point.
These expansion scopes are far below total wall time, so attributing the
remaining iron runtime primarily to these coupled HVP expansions would be
incorrect without further stage evidence.

Archived four-step job 1248959 gives singlet 4/50, 529.15 s,
-2244.7182260015 Eh and triplet 4/41, 610.02 s,
-2244.7615530749 Eh. Several numerical-policy revisions separate these
runs, and nodes differ. They are historical trajectory comparisons, not an
isolated timing test of metric factorization. Four-step energies cannot
establish endpoint agreement with converged PySCF results.

Raw traces and GNU-time records remain in
`fecl2-ritz-cost-1249408/native/` and `coarse-regression-1249412/`.
`summarize_ritz_cost.py RUN_DIRECTORY` extracts the scopes/action counts.
