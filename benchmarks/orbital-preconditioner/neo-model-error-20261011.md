# NEO model error: separating local orbital and structure curvature

## Correction to the previous attribution

The admission patch is a single backend-selection change, but it changes more
than one mathematical contribution. The core-only branch disables both
`local_active_response` and `structure_response`. The coupled branch calls
`apply_orbital_coupling`, which includes local-active orbital curvature,
without internally solving structure response; it supplies structure forcing
to a separate coupled solver. Thus the earlier admission experiment isolates
the model switch, not structure coupling alone. Previous reports have been
corrected to state this explicitly.

## Operator decomposition

**Correction from the subsequent same-point audit:** local-active alone is
not the fixed-structure orbital Hessian. Coupled orbital action additionally
includes the selected-state overlap-normalization gauge response. The earlier
formula omitted this contribution; the corrected decomposition follows.

Let $A_0$ denote direct-core and fixed-upstream pullback curvature, $L$ the
computational local-active contribution, $G$ the overlap-normalization gauge
contribution, $B$ the horizontal orbital/structure coupling, and $C$ the
structure tangent operator. At a stationary structure solution,
with the inverse restricted to the appropriate tangent space, the formal
relaxed orbital operator is

$$
\bar A=A_0+L+G,\qquad H_{\mathrm{rel}}=\bar A-B^T C^\dagger B.
$$

The production core-only approximation is $A_0$, not the complete
fixed-structure orbital operator $\bar A$. The omitted response is therefore
the combination

$$
\Delta H=L+G-B^T C^\dagger B,
$$

not solely the structure Schur term. If a core-only shifted solve supplies
a candidate $p$, its full relaxed residual includes $\Delta Hp$. A separate
small core-only KKT residual cannot bound this omitted term. Singular or
degenerate structure response requires the operator's tangent projection
and compatibility conditions; this notation does not justify an unrestricted
inverse of $C$.

The new experiment leaves structure coupling disabled and adds only $L$.
It is a controlled toggle in the same orbital-only solver. Comparison with
coupled NEO also changes solver representation/workspace; trajectories alone
must not be read as a direct measurement of the Schur term.

## Existing accepted-trajectory evidence

NEGPIA core-only first shows rejected trials at point 4: two rejects, an
accepted radius of 0.0625 rather than 1, and 40 orbital HVPs for that point.
Further poor ratios and rejections leave its accepted radius at 0.0009765625
by point 14. The coupled trajectory has no rejected trials and keeps radius
1 throughout all six points. Its final three actual/predicted decrease ratios
are approximately 0.962, 0.996 and 1.000.

Other core-only systems generally retain radius 1 and do not reject trials.
Their weaker late-stage gradient reduction can add points, but there is no
comparable radius collapse. Ratios near one do not prove a correct full
gradient model: they measure a scalar along an accepted direction only.
Physical gradient norms at different charts are not gradient difference
vectors. No transported gradient-error claim is made from these logs.

## Local-active-only ablation

`neo-local-active.patch` changes only the local-active flag in the core-only
branch. Production tolerances, structure admission, geometry and structure
response remain unchanged. Build and NEGPIA runs use Hanhai25 job `281633`;
the remaining five matched pairs use parallel array `281634`. Each pair uses
the same node, historical input, 32 compact-affinity threads, exact Libcint
and Davidson. GNU time supplies native peak process RSS. The source mirror
and production executable are restored after the experimental build.

Each cell is **current core-only / local-active-only**:

| System | Points | Orbital HVP | Rejected trials | Total s | RSS MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| 240 | 8/20 | 12/77 | 0/5 | 2.166/6.559 | 1688.5/1733.7 |
| 67970 | 7/14 | 13/67 | 0/4 | 3.226/7.977 | 1977.0/1989.7 |
| CONYAH | 7/18 | 13/84 | 0/7 | 2.001/6.289 | 1617.4/1669.5 |
| QOCGUM | 7/14 | 15/106 | 0/7 | 9.136/28.827 | 2083.5/2128.9 |
| YEPKAF | 7/31 | 25/273 | 0/7 | 9.247/44.614 | 2299.5/2273.5 |

All ten panel runs explicitly converge and exit zero. Their final energy
differences are below 4.2e-8 Eh for the first four cases; YEPKAF differs by
2.64e-7 Eh. All use the same outer criteria, not exactly identical endpoint
accuracy. Native run logs and hashes are retained. No profiler overhead is
included and no precise hotspot attribution is inferred from wall time.

NEGPIA local-only is still running. The saved CSV is an explicitly incomplete
snapshot at 106 points and 1126 HVPs, with no final time, RSS or energy filled
in. A later scheduler read reached point 124. It has not been declared
converged or failed, and no new optimizer iteration cap was imposed.

## What the evidence supports

Adding local-active curvature alone is not a viable improvement on these
cases. The same-point follow-up distinguishes the overlap-normalization
gauge term from horizontal structure relaxation. The apparent large
cancellation in the computational local/structure split must not be attributed
solely to the physical Schur term. See `neo-component-audit-20261011.md`.

Before implementing a model-error-driven admission policy, measure the two
contributions on exactly the same accepted point and actual candidate direction:
local-active image, structure-Schur image, and their sum. Validate their sum
against a same-chart relaxed-gradient directional difference. This must retain
normalization, tangent projection and gauge consistency, and should be a
bounded diagnostic on NEGPIA rather than another full multi-system optimization.
Such measurements distinguish a real omitted-response error from a derivative
inconsistency or a subproblem/globalization difference.

No production strategy changes or new empirical thresholds are introduced.
Data are in `results/history/neo_local_active_20261011.csv`; raw logs are under
`results/studies/locality-20261011/neo-admission-281633/` and
`results/studies/locality-20261011/neo-admission-panel-281634/`, with corresponding
remote directories in `/home/guqqgroup/taoxia/xmvb-runs/native-20261011/`.
