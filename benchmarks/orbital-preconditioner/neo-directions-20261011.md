# NEO orbital expansion and preconditioner study

Reference: `fbb7fc26`, frozen executable SHA256
`68a315316d83f3b3923f9dd3ed9b94fb08588669078d763de8316c73d416c53b`.
All workloads and builds run on Hanhai25 Slurm compute nodes, with 32
compact-affinity OpenMP cores, BLAS threads fixed to one, Davidson,
exact Libcint integrals, gradient tolerance 1e-3 and energy tolerance 1e-7.
No optimizer work cap or tolerance relaxation is introduced.

## Historical strategy and the mathematical difference

The historical and current one-electron preconditioners use the same
positive spectral completion of each local normalized-orbital curvature:

$$
P_i=V_i\operatorname{diag}
\left[\max\left(|\lambda_j|,\sqrt{\epsilon}\max_l|\lambda_l|\right)\right]V_i^T.
$$

Here $V_i$ and $\lambda_j$ are the local surrogate eigenvectors and
eigenvalues, and $\epsilon$ is machine precision. This is a surrogate for
preconditioning, not the exact orbital Hessian. The exact projected model
retains negative curvature. The regularization was already present in
`c2776c4`; it is not a newly tuned parameter explaining the regression.

The historical coordinate residual test does not enforce today's physical
KKT certificate. FeCl2 spends 60 HVPs on the first step in both versions,
but its second step increases from five to 226 HVPs. Previously measured
true nonlinear gradients show that short second-point trials do not resolve
physical stationarity. Restoring the old stop alone is not an equal-accuracy fix.

NEGPIA's 45-to-51 HVP increase is a trajectory change, primarily in the
middle of the optimization, not excess final-response certification. The
restored coupled model and conservative forcing yield six accepted points
instead of seven. Matching historical iteration counts is not sufficient
to establish wall-time improvement.

## Rejected orbital expansions

Let $Q$ be the retained Euclidean-orthonormal orbital basis and $r$ the
orbital correction right-hand side. Ordinary Davidson expands by
orthogonalizing $P^{-1}r$ against $Q$.

The first experiment instead solves the constrained positive model
$Pd+Q\mu=r$, $Q^Td=0$:

$$
d=P^{-1}r-P^{-1}Q(Q^TP^{-1}Q)^{-1}Q^TP^{-1}r.
$$

The second uses accepted-point exact orbital images $AQ$, not transported
secant history, to complete the inverse surrogate:

$$
D=Q-P^{-1}AQ,\qquad
H_P=P^{-1}+D(D^TAQ)^\dagger D^T.
$$

Neither introduces extra HVPs for its algebraic update. Both nevertheless
generate worse orbital subspaces on the target systems. The native evidence
is decisive; neither implementation, auxiliary cache nor interface remains.

| Expansion | Jobs | FeCl2 HVP: reference / candidate | FeCl2 wall s | NEGPIA HVP | NEGPIA wall s |
| --- | --- | ---: | ---: | ---: | ---: |
| Constrained inverse | 281794 / 281797 | 287 / 494 | 6.098 / 17.024 | 51 / 63 | 9.831 / 11.264 |
| Inverse SR1 completion | 281807 / 281808 | 287 / 338 | 5.999 / 10.166 | 51 / 52 | 9.654 / 9.785 |

Pairs share one node and identical input hashes. The constrained inverse
also fails the MnF2 response certificate after 23 accepted points. The SR1
completion passes the unit tests and six operator checks, but its measured
target performance does not justify retention.

## Analytic diagonal: preserve the formula, remove unused transforms

The analytic model previously stored $(pq|ab)$ and $(pa|qb)$ for all
$p,q$, where $a,b$ denote active orbitals. Inspecting the diagonal formula
shows that its external contractions need only

$$
J_{pab}=(pp|ab),\qquad K_{pab}=(pa|pb),
$$

plus the active tensor $(ab|cd)$. The generalized orbital-gradient term
needs only its diagonal, not every $g_{pq}$. Its active diagonal is

$$
g_{aa}=\sum_u(h_{au}+V^{\rm core}_{au})\gamma_{ua}
       +\sum_{uwz}(au|wz)\Gamma_{wzua}.
$$

Here $h$ is the one-electron Hamiltonian in the accepted orthogonal
computational frame, $V^{\rm core}$ its inactive Coulomb/exchange potential,
$\gamma$ the active one-particle density, and $\Gamma$ the spin-free
two-particle density used by the original analytic model. This frame is
only a way to construct the preconditioner; nonorthogonal VB orbital
updates and exact HVPs are unchanged.

The two external tensors therefore shrink from $2n_{\rm bf}^2n_a^2$
entries to $2n_{\rm bf}n_a^2+n_a^4$. External pair contraction work drops
from $O(N_{\rm AO-pair}n_{\rm bf}^2n_a^2)$ to
$O(N_{\rm AO-pair}n_{\rm bf}n_a^2)$, excluding graph application and
active-only contractions. AO graph work to produce the mixed pair images
is unchanged. The exact and RI representations share this slicing logic.

Workers share read-only AO pair blocks and own only active-pair vectors or
active-square contraction slices. No new per-thread integral or large-matrix
replica is introduced. Existing AO-kernel workspace behavior is outside
this claim; the change does not prove that all memory replication is removed.

Job 281820 passes nine unit executables and six F2 HAO/OEO/state-averaged
operator checks in exact/RI modes. The analytic diagonal finite-difference
error is 5.99e-8 with the original finite-difference acceptance tolerance.
Default one-electron trajectories are unchanged. Full-diagonal ablations
are jobs 281796 (original) and 281822 (compact).

Additional explicit analytic-diagonal optimization checks are recorded in
job 281985: original/retained pairs for F2 HAO, OEO and two-state averaging
in exact and RI modes. All 12 converge, with matching within-pair initial
energies and final-energy differences below 4.1e-9 Eh. These are numerical
smoke checks, not native speed measurements. The earlier exploratory job
281975 must also be retained: its F2 OEO exact attempt exhausted the trust
radius and stopped the old fail-fast script. The paired rerun does not erase
that failure or establish deterministic/global convergence of the optional
diagonal preconditioner. The default remains one-electron.

| System | Original diagonal | Compact diagonal | Interpretation |
| --- | ---: | ---: | --- |
| NEGPIA | 6 points, 24 HVP, 248.993 s | 6 points, 24 HVP, 164.369 s | Both converge; remaining transform cost exceeds the 9.8–10.0 s one-electron run. |
| FeCl2 | 91 points, 3022 HVP, 522.375 s | 98 points, 2906 HVP, 184.942 s | Both fail the requested structure-response certificate; not converged benchmarks. |
| MnF2 | 25 points, 573 HVP, 82.579 s | 24 points, 572 HVP, 21.478 s | Both fail to converge; not speedup evidence for an accepted optimizer. |

Each representation comparison reused the same respective case node,
but allocations were separate. Arithmetic reassociation can change raw
orbital gauges, coordinate infinity norms and failed nonlinear trajectories;
these rows are not claims of identical failed paths or universal optimizer
dominance. NEGPIA's converged endpoint energies agree within 2e-12 Eh and
physical gradient norms within 1e-12.

## Rejected active-only analytic preconditioner

An additional, system-independent sector decomposition kept the original
one-electron blocks for inactive orbitals and used the analytic model for
active orbitals. In the latter chart, physical directions are first projected
out of the inactive span. The inactive--external exchange slices therefore
do not enter the active rotation submatrix. A unit test confirmed that
zeroing these slices leaves the active submatrix unchanged to 1e-13.

This removes a real construction cost, but it is not a successful orbital
preconditioner. Job 281835 passes all correctness checks; array 281836
runs one-electron / active-only pairs on identical nodes and inputs:

| System | Points | HVP | Wall s | Peak RSS MiB |
| --- | ---: | ---: | ---: | ---: |
| FeCl2 | 3 / 3 | 287 / 364 | 5.956 / 8.198 | 389.1 / 412.9 |
| NEGPIA | 6 / 6 | 51 / 37 | 9.718 / 17.639 | 3776.6 / 3830.6 |
| MnF2 | 14 / 21 | 181 / 341 | 3.894 / 9.969 | 496.8 / 506.4 |
| 67970 | 5 / 5 | 11 / 13 | 3.334 / 7.051 | 1960.0 / 1986.8 |
| QOCGUM | 5 / 5 | 13 / 15 | 6.814 / 11.257 | 2076.7 / 2146.0 |
| 240 | 5 / 5 | 10 / 14 | 1.930 / 4.404 | 1691.5 / 1741.8 |

All converge, but none establishes a wall-time improvement. NEGPIA's
fewer HVPs do not repay the construction. The enum, CLI option, partial
factory mode and chart branch have all been deleted. Default selection
remains one-electron, with no molecule-specific preconditioner rule.

These results do not establish that the existing preconditioner is erroneous.
They do establish that changing only projection algebra, completing only
measured orbital images, or replacing local active blocks with the rotation
diagonal is not a demonstrated repair of the target regressions.

## Retained candidate-local Schur seed

Unconditionally initializing every strict solve with $Cz=-BQ_0$ reduces
FeCl2 to 235 HVPs, but wastes structure actions elsewhere: QOCGUM rises
from 72 to 285 structure actions without reducing its 13 HVPs. Jobs
281948/281949 record this experiment; the unconditional policy is removed.

The retained policy uses the existing readiness and certification residuals:

$$
\tau_c=\epsilon_{\rm abs}+\frac{\|g\|_2^2}{1+\|g\|_2},\qquad
\tau_p=\epsilon_{\rm abs}+\eta\|g\|_{G^{-1}}.
$$

Here $G$ is the retraction metric, $g$ the accepted-point orbital gradient,
and $\eta$ the existing forcing term. In a strict solve with an empty
structure basis, the seed is admitted only when

$$
\sqrt{\|r_o\|_2^2+\|r_s\|_2^2}\leq\tau_c,
\qquad
\sqrt{\|r_o\|_{G^{-1}}^2+\|r_s\|_2^2}>\tau_p.
$$

These conditions identify a coordinate-ready but physically unresolved
candidate. They are a basis-expansion trigger, **not** a coordinate
convergence certificate. No new empirical tolerance is introduced.

The scalar/phase trace in job 281969 shows the concrete change. At FeCl2's
second point, after five HVPs the empty-response model reports physical
orbital residual 2.2939e-3 and structure residual 2.4079e-4. Retaining the
candidate response adds no HVP and exposes orbital residual 4.1677e-3,
while reducing the structure residual to 6.6653e-6. The larger corrected
orbital residual is useful information: subsequent directions refine the
coupled model instead of continuing under an underestimated response.
The second-point orbital work drops from 226 to 185 HVPs; the first and
third points remain 60 and one. The final-backreaction checks at the second
point drop from 18 to 13, in addition to the one seed.

Instrumented scopes are inclusive, not additive. In job 281969 the NEO
solve takes 3.310 s, orbital HVPs 2.495 s, projected solves 0.498 s and
the 15 structure solves 0.048 s. These are phase measurements from the
instrumented run, not native speedup timings. Orbital HVP work remains
the main measured cost; optimizing only the structure-solver kernel would
not address it.

With the current structure basis empty, the cached candidate has $r_s=Bp$.
Compute one response $q_0=-C^\dagger Bp$ using the existing residual-certified
structure solver, and retain its coupled image $(B^Tq_0,Cq_0)$ as a basis
direction. All finite-precision response tolerances are unchanged. After
reprojection, the structure coefficient varies with the new candidate;
one retained response does not claim to solve all future $Cq=-Bp$ equations.
The original complete coupled residual and final response-backreaction check
remain mandatory. An accepted-point flag prevents repeated zero/degenerate
seed attempts. Ordinary early nonlinear trials remain deferred.

No new HVP, complete coupling matrix, integral copy or pair cache is required
for the seed. It adds one existing vector/image column, of size
$O(n_o+n_s)$, only when the mismatch requires it. The exact orbital
preconditioner and Hessian actions are not modified.

Jobs 281956 and 281958 build/check and run the final six-case panel. The
new anisotropic-metric unit test is rebuilt/run in job 281967. It requires
the seed before the second orbital image and compares the converged step
with the explicit coupled Newton solution. The original small-response,
early-trial deferral, strict residual and bounded-subspace tests remain
unchanged. Nine unit executables and six exact/RI operator checks pass.

| System | Points | HVP | Structure actions | Wall s | Peak RSS MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| FeCl2 | 3 / 3 | 287 / 246 | 1008 / 784 | 5.933 / 5.154 | 399.5 / 391.4 |
| NEGPIA | 6 / 6 | 51 / 51 | 70 / 69 | 9.613 / 9.639 | 3774.2 / 3739.7 |
| MnF2 | 14 / 14 | 181 / 181 | 426 / 424 | 4.164 / 4.132 | 491.6 / 477.1 |
| 67970 | 5 / 5 | 11 / 11 | 61 / 61 | 3.264 / 3.271 | 2000.3 / 2056.0 |
| QOCGUM | 5 / 5 | 13 / 13 | 72 / 72 | 7.039 / 7.001 | 2057.8 / 2118.3 |
| 240 | 5 / 5 | 10 / 10 | 43 / 43 | 1.943 / 1.860 | 1741.2 / 1725.9 |

Rows show reference / retained candidate; each pair shares a node, input
hash and settings. The candidate SHA256 is
`95058dd2331f7d0250989136aa8058d598899938e980d77cd4ed506a93dd0b8d`.
The FeCl2 improvement is 14.3% fewer HVPs and 13.1% lower total wall time,
not a return to the historical 72 actions. NEGPIA's regression is unresolved;
unchanged action counts and small timing fluctuations are not a speedup.

Accuracy criteria are unchanged, but endpoints must not be described as
identical. FeCl2 changes from -2181.617645215964 to -2181.617644104857 Eh,
an increase of 1.1111e-6 Eh. Its physical gradient changes from 9.8857e-4
to 9.4766e-4; its final accepted energy change is 3.0273e-9 Eh. Both pass
the same stopping tests. Equal residual tolerances do not prove equal
energy errors in a soft-curvature region; this endpoint difference remains
part of the unresolved convergence analysis. No tolerance was widened to
obtain the speed improvement.

The retained changes do not establish universal dominance or a complete
solution of the physical-metric orbital refinement problem. Raw native
evidence, including all failed runs and prototype hashes, is recorded in
`results/history/neo_directions_20261011.csv`.
