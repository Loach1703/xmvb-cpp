# Accepted-point normalization gauge

## Exact factorization

Let $V$ contain the selected structure coefficients, $E$ their diagonal
energy matrix, and $w_j$ the equal positive state weights. Let $H$ and $M$
be the accepted structure Hamiltonian and overlap operators. Retain the
accepted Ritz residual $R=HV-MVE$; it is not set to zero. Define

$$
G=V^T M V,\qquad D(p)=V^T\delta M[p]V,
\qquad Q(p)=V^T\bigl(\delta H[p]V-\delta M[p]VE\bigr).
$$

The normalization gauge and its selected-space multiplier are

$$
\Gamma(p)=-\frac12 G^{-1}D(p),
$$

$$
\Lambda(p)=G^{-1}\left[
Q(p)+G\bigl(E\Gamma(p)-\Gamma(p)E\bigr)
+V^TR\Gamma(p)\right].
$$

These are the same equations used by the coupled action. $G^{-1}$ denotes
a small LU solve, not a structure-overlap inverse. The selected metric
normalization check remains unchanged.

Let $K(\delta V,L)$ denote the existing structure-response adjoint pulled
back to reduced orbital coordinates. For the matrix unit $e_{ij}$, cache

$$
a_{ij}=K(Ve_{ij},0),\qquad m_{ij}=K(0,e_{ij}).
$$

Adjoint reciprocity gives

$$
Q_{ij}(p)=\frac{a_{ij}^Tp}{2w_j},\qquad
D_{ij}(p)=\frac{m_{ij}^Tp}{w_j}.
$$

The gauge action can therefore be evaluated as

$$
\mathcal Gp=\sum_{ij}\left[
\Gamma_{ij}(p)a_{ij}-\Lambda_{ij}(p)m_{ij}\right].
$$

Finally, with $A_0$ the direct-core/fixed-upstream contribution and $L_o$
the computational local-active contribution,

$$
\bar A p=(A_0+L_o)p+\mathcal Gp.
$$

This is a representation change, not a truncation of the gauge. It does not
solve or approximate the horizontal structure response. Coupled NEO retains
its existing $B$, $C$, residual checks and response admission policy.

## Cached-channel experiment: not retained

At each accepted point, initialize two orbital adjoint channels per ordered
selected-state pair. Subsequent gauge actions use dot products and linear
combinations of those channels, plus small selected-state matrix solves.
Storage and per-action vector work are $O(n_p n_s^2)$, where $n_p$ is the
reduced orbital dimension and $n_s$ the number of selected states. There is
no unique-spin-pair-square cache, orbital Hessian, or per-thread integral
copy. The channels belong to one operator/accepted-point/chart workspace;
they are not reused across changed orbital points.

The initialization cost is $2n_s^2$ adjoint applications. Thus this removes
repeated work, not all accepted-point work; it need not win when very few
directions are applied or many states are selected. The full local-active
orbital HVP and integral transformations still occur for each direction.

This exact representation was implemented and certified against two
independent directions for HAO, OEO and two-state-average F2, in both exact
and RI modes. Errors were below $6\times10^{-16}$. Nevertheless, all six
native cases ran slightly slower. Jobs 281665/281666 measured the first
version; 281672/281673 reused the single-state accepted gradient;
281674/281679 additionally removed zero-response CI operations. Initialization
and loss of the existing fused backward outweighed the eliminated forward
work on this panel. The cached implementation and its storage were deleted.
It is not a production alternative or fallback.

## Retained implementation

`ExactHvpOperator::apply_orbital` now uses the existing fused local/gauge
backward. It requests only the gauge from
`prepare_equal_weight_generalized_eigen_coupling`, obtaining $D$ and $Q$
from the already computed directional selected-state images. It omits the
full horizontal forcing, its gauge-shifted image, selected-span projection
and scaled RHS. The accepted Ritz residual and normalization checks remain.

The non-coupled NEO branch uses this complete orbital operator instead of
the incomplete computational core-only model. The coupled branch remains
unchanged. There is no new accepted-point channel cache or dense Hessian.

The direct-CI response adjoint also skips exterior transforms, sigma and
integral adjoints for **exactly zero** coefficient responses. Overlap
multiplier contributions are still evaluated. States are consumed one at a
time instead of retaining all their transformed response carriers together.
No numerical threshold or approximate truncation is introduced.

## Validation

`check_exact_ctx_hvp` compares two independent directions of the orbital-only
operator with the orbital image of `apply_orbital_coupling`. Its normalized
maximum-element difference must be at most $10^{-10}$. The reference still
constructs directional selected-state images and the full gauge adjoint.
The retained path forms only selected-space projections before the same fused
backward. A separate synthetic unit test uses two normalized roots with a
deliberately nonzero Ritz residual and checks that both full forcing fields
remain empty in gauge-only mode.

Hanhai25 job 281684 builds and checks HAO, OEO and two-state-average F2 for
both exact and RI integrals, then compares NEGPIA with the preceding complete
orbital-gauge executable. Array 281685 compares five independent historical
systems after the checks finish. Job 281686 adds the synthetic unit test and
captures opt-in phase timing separately; array 281691 adds MnF2 and FeCl2.
Job 281693 rebuilds the final source after renaming the former `core` result
to `orbital` and repeats all eight unit executables and six operator checks.
All pass; the maximum orbital-only comparison error is $7.08\times10^{-16}$.
Its final executable SHA256 is
`42e74d44245d68130f58ac592a65087e97bcc64a701d1dbccf34d58a351a8d03`.
Its native NEGPIA pair is 8.820/8.768 s, again with 7 points and 45 HVPs.
Each native pair uses the same node, input,
32 compact-affinity OpenMP threads and one BLAS thread. Convergence thresholds
remain $10^{-3}$ for the projected gradient and $10^{-7}$ Eh for adjacent
energies. Profiling runs are not pooled with native timing samples.

## Native six-case results

Each cell is previous complete orbital-gauge / retained fused gauge-only.
These are single matched native samples, not statistical speed estimates.

| System | Points | HVP | Total s | RSS MiB |
| --- | ---: | ---: | ---: | ---: |
| NEGPIA | 7/7 | 45/45 | 8.711/8.778 | 3798.0/3813.1 |
| 240 | 8/8 | 13/13 | 2.300/2.304 | 1697.5/1682.1 |
| 67970 | 7/7 | 13/13 | 3.653/3.652 | 2037.7/2025.6 |
| CONYAH | 7/7 | 13/13 | 1.999/2.045 | 1719.6/1661.0 |
| QOCGUM | 7/7 | 14/14 | 7.584/7.518 | 2100.5/2094.9 |
| YEPKAF | 7/7 | 23/23 | 9.346/9.207 | 2250.6/2218.4 |

All 12 applications converge with no rejected trials and zero horizontal
response actions. Endpoint energies agree within $2\times10^{-12}$ Eh.
The operator checks differ by at most $4.72\times10^{-16}$.
This removes unnecessary RHS construction without a meaningful end-to-end
speedup or memory reduction on these small cases. The improvement over the
former incomplete core model is the already established NEGPIA trajectory
recovery (17 to 7 points, 168 to 45 HVPs), not a new timing gain from caching.

## Phase diagnosis

The additional matched native transition-metal checks also converge:

| System | Points | HVP | Structure actions | Total s | RSS MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| MnF2 | 14/14 | 180/180 | 437/437 | 4.164/4.139 | 481.7/494.0 |
| FeCl2 | 3/3 | 287/287 | 1006/1007 | 6.180/6.235 | 402.1/399.5 |

Their energies agree within $5\times10^{-12}$ Eh. These exercise the existing
coupled branch and the exact-zero response-adjoint change, not just the
orbital-only branch. The one-action FeCl2 response difference can arise from
iterative arithmetic near a stopping boundary; no solver criterion was changed.
The selected inputs are saved historical inputs, whose initial guesses differ
from other transition-metal benchmark series; counts and times should not be
pooled across those series.

Job 281686 uses opt-in phase scopes in a separate matched NEGPIA pair.
Scopes can nest, so inclusive totals must not be summed. The 45 HVPs take
4.299/4.300 s (previous/retained). Selected-direction preparation takes only
0.070/0.069 s, explaining why eliminating horizontal RHS overhead cannot
produce a large native gain here. In the retained run, active integral
direction preparation takes 1.985 s, AO response 1.335 s and orbital pullback
0.780 s. These dominate selected-direction work. Future substantial runtime
optimization should investigate those integral transformations and AO
contractions, not add another normalization-gauge cache.

Eight unit executables pass, including the new finite-Ritz gauge test
(error $4.39\times10^{-17}$ at Ritz residual norm $3.14\times10^{-3}$).
The final source operator checks also pass for all six orbital/integral
combinations. Native six-case records, input hashes and executable hashes
are in `results/history/neo_orbital_fused_20261011.csv`; raw logs are under
`results/studies/locality-20261011/` and the Hanhai25 native run directory.
