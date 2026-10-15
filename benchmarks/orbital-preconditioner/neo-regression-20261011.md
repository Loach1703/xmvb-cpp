# NEO regression against c2776c4

The iteration/action regression is repaired for 67970 and QOCGUM: both again
take five accepted orbital points and respectively 11 and 13 orbital HVPs.
FeCl2 still takes 287 HVPs. The experiments below do **not** establish a safe
way to recover its historical 72-HVP result at the current physical accuracy.

## Evidence

The historical reference is `c2776c4`, including its raw `optimizer.tsv` and
the corrected accepted-point counts in `results/history/neo_optimized_20261011.csv`.
It was not rebuilt or rerun in this investigation. Historical wall times are
therefore context, not matched speed measurements.

The fresh reference is `ac9630f2`, executable SHA256
`42e74d44245d68130f58ac592a65087e97bcc64a701d1dbccf34d58a351a8d03`.
The retained candidate has SHA256
`68a315316d83f3b3923f9dd3ed9b94fb08588669078d763de8316c73d416c53b`.
Job 281771 builds and checks the candidate and runs the NEGPIA pair.
Array 281772 runs the other five pairs concurrently. Every pair shares a node,
input hash, 32 compact-affinity OpenMP cores, and one BLAS thread. All use
Davidson and exact Libcint integrals. Build flags are
`-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`. No optimizer work cap was added.

All 12 native applications converge. Eight unit executables pass, as do
the two-direction orbital/coupled HVP comparisons for F2 HAO, OEO and
two-state averaging in both exact and RI modes. These are derivative/operator
checks; no RI optimization speedup is claimed.

Machine-readable results, input/binary hashes and remote raw-log locations:
`results/history/neo_regression_20261011.csv`.

## Two retained changes

### Action work, not obsolete cofactor work

The old revision always admitted a nonempty structure tangent. The fresh
reference rejected it for 67970 and QOCGUM. Its work estimate counted only
AO-pair rows for the exact core, while estimating structure response using
the former cofactor sandwich even when the actual action was direct CI.

Let $P_a=n_a(n_a+1)/2$ be the number of packed active-orbital pairs,
$E_{\rm AO}$ the number of retained AO-pair graph entries, and $U_\alpha,U_\beta$
the unique-spin-string counts. The relevant leading work is

$$
W_{\rm core}\sim E_{\rm AO}P_a,\qquad
W_{\rm direct\,CI}\sim n_s U_\alpha U_\beta P_a^2.
$$

Here $n_s$ is the selected-state count. The core estimate now uses the actual
CSR entry count, or the leading two-direction count for compact integrals.
The structure estimate follows the retained action representation. These
are leading work estimates, not complete wall-time predictions: response
iterations, exterior transforms and parallel efficiency also matter.
The RI graph-empty core proxy is unchanged.

No graph scan, Hessian, integral copy or unique-pair cache is added.

### Preserve the historical forcing budget

The first accepted points agree with the historical trajectories. The next
point diverged because the larger physical gradient norm made the absolute
inner budget wider, admitting a one-direction model where the old strategy
continued to a second curvature direction.

Write $g_c=\|g\|_2$ for the coordinate gradient norm and
$g_p=\sqrt{g^TG^{-1}g}$ for its physical dual norm. $G$ is the existing
retraction metric, not a newly assembled matrix. The existing forcing rule is

$$
\eta(t)=\max\left(\epsilon_{\rm round},\frac{t}{1+t}\right).
$$

For nonlinear trial candidates, the relative physical forcing now obeys

$$
\eta_{\rm trial}=
\min\left[\eta(g_p),\frac{\eta(g_c)g_c}{g_p}\right],
\qquad
\tau=\epsilon_{\rm abs}+\eta_{\rm trial}g_p.
$$

The trial-enabled branch already requires a positive physical gradient,
so there is no division at zero. Thus the absolute budget cannot exceed
either the historical coordinate budget or the current physical budget.
No empirical coefficient is introduced. The existing near-stationarity
branch continues to use the physical budget alone. The final gradient,
energy, actual-trial and structure-response certificates are unchanged.
Trial readiness remains distinct from certified KKT convergence.

## Native results

Cells show fresh reference / retained candidate. Points count every accepted
orbital update, not macro groups. Times are end-to-end, including setup/output.

| System | Points | Orbital HVP | Structure actions | Wall s | Peak RSS MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| 67970 | 7 / 5 | 13 / 11 | 0 / 61 | 3.497 / 3.327 | 1981.8 / 2052.2 |
| QOCGUM | 7 / 5 | 14 / 13 | 0 / 72 | 8.096 / 7.393 | 2043.2 / 1999.2 |
| FeCl2 | 3 / 3 | 287 / 287 | 1008 / 1009 | 6.216 / 6.188 | 401.4 / 400.3 |
| 240 | 8 / 5 | 13 / 10 | 0 / 43 | 2.188 / 1.906 | 1786.9 / 1782.3 |
| MnF2 | 14 / 14 | 180 / 181 | 437 / 425 | 3.954 / 4.048 | 488.0 / 502.8 |
| NEGPIA | 7 / 6 | 45 / 51 | 0 / 69 | 8.952 / 9.814 | 3813.3 / 3778.1 |

This is **not** universal dominance: NEGPIA takes fewer points but more HVPs
and longer wall time; MnF2 is slightly slower. The historical comparison is:

| System | Historical points / HVP | Historical wall s | Candidate points / HVP | Candidate wall s |
| --- | ---: | ---: | ---: | ---: |
| 67970 | 5 / 11 | 3.165 | 5 / 11 | 3.327 |
| QOCGUM | 5 / 13 | 6.330 | 5 / 13 | 7.393 |
| FeCl2 | 4 / 72 | 5.285 | 3 / 287 | 6.188 |
| 240 | 5 / 10 | 10.681 | 5 / 10 | 1.906 |
| MnF2 | 14 / 199 | 9.860 | 14 / 181 | 4.048 |
| NEGPIA | 6 / 51 | 13.416 | 6 / 51 | 9.814 |

QOCGUM's candidate endpoint agrees with the historical energy within
$4\times10^{-12}$ Eh; 67970 within $6\times10^{-10}$ Eh. Their historical
action/point trajectories are recovered, but historical wall-time parity is
not yet demonstrated. Different allocations and added certificates must not
be hidden behind matching iteration counts.

## Why FeCl2 went from 72 to 287

Both revisions spend 60 HVPs on the first step and produce the same energy
decrease to about $2\times10^{-11}$ Eh. The difference is the second step:
5 HVPs historically, 226 in the fresh reference. The historical solve tested
coordinate residuals; the fresh solve must also resolve physical stationarity
and the actual orbital backreaction of its structure residual.

The historical endpoint is $-2181.617639896249$ Eh, with a reported physical
gradient norm $1.1490\times10^{-3}$. The candidate endpoint is
$-2181.617645215969$ Eh, with physical norm $9.8859\times10^{-4}$.
The old physical-metric implementation also differs; its reported norm is
not a remeasurement of its final orbitals under today's metric. Restoring
the old stopping behavior is therefore not an equal-accuracy repair.

Job 281748 audits real, fully relaxed nonlinear trials without accepting them:

| Second-point orbital HVPs | Model orbital residual | Actual trial physical gradient |
| ---: | ---: | ---: |
| 5 | 0.002294 | 0.003863 |
| 50 | 0.001314 | 0.004098 |
| 100 | 0.001696 | 0.003753 |

All three trial energies decrease, but none satisfies physical stationarity.
Energy decrease alone cannot certify these short solves. The scalar trace in
job 281715 also shows 18 structure-response certificates at the second point;
the first changes the residual from 0.000997 to 0.004152.

The instrumented scopes in that job are inclusive and must not be summed:

| Scope | Calls | Seconds |
| --- | ---: | ---: |
| NEO solve | 3 | 4.301 |
| Orbital HVP | 287 | 3.029 |
| Projected NEO solve | 304 | 0.824 |
| Structure solve | 19 | 0.063 |

The structure solver itself accounts for about 1.5% of the NEO solve time,
not the main runtime bottleneck despite its large action count. The evidence
points to inefficient orbital-space refinement under the physical certificate;
it does not identify a broken integral contraction or justify loosening accuracy.

## Rejected implementations

All of these were built and tested on Hanhai25; none remains in production:

| Experiment | Native job | FeCl2 points / HVP | Matched wall s: reference / experiment |
| --- | ---: | ---: | ---: |
| Enable cheap tail trials with any true-gradient progress | 281697 | 56 / 372 | 5.536 / 25.448 |
| Replace tail directions with metric Riesz directions | 281704 | 3 / 412 | 5.446 / 10.026 |
| Add Riesz directions alongside block-preconditioned directions | 281712 | 3 / 326 | 6.141 / 7.294 |
| Old expansion priority plus paired response correction | 281726 | 3 / 325 | 6.176 / 6.936 |
| Close candidate response earlier | 281736 | 3 / 279 | 6.150 / 6.742 |
| Physical-norm recycling epochs | 281742 | 3 / 347 | 5.563 / 6.478 |
| Cached retained-response Woodbury preconditioner | 281762 | 3 / 322 | 6.100 / 7.124 |

These experiments change algebraic strategy, not fitted tolerances. They show
why neither blind restoration of cheap trials nor a structure-only correction
is a demonstrated solution to FeCl2's orbital refinement cost. The next study
must quantify the missing orbital curvature in the existing preconditioner,
with the full physical residual and actual nonlinear gradient as references.
