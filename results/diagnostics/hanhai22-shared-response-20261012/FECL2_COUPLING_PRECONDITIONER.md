# FeCl2: curvature cancellation and sampled-coupling preconditioning

## Evidence, not a convergence claim

Restored-path diagnostic 1248906 converged in 16 macros / 334 HVPs,
11.90 s GNU wall, 526.9 MiB RSS, at -2181.624464775997 Eh.
All five unit tests passed. The trace separately pairs orbital, structure
and coupling images already computed by each coupled HVP; no extra HVP is
needed. Forward/backward coupling pairings agree to at most 2.193e-15 in
the tail, after scaling by max(1, their magnitudes).

For one admitted direction in macro 12, the orbital quadratic contribution
is 0.1272933610, the structure contribution 0.1120552579, and the sum of
both coupling contributions -0.2249607196. The total retains only 6.01%
of the uncoupled curvature. In macro 15 another direction retains 7.85%.
Thus some sampled soft directions arise from cancellation between blocks,
which an independent orbital/structure diagonal inverse does not represent.
This does not prove that all slow modes have this origin.

At nonlinear keyframes, 90.28%--99.62% of the squared physical-dual gradient
correction in macros 9--16 lies in the orbital coordinates. The gradient
correction is not primarily a structure-norm accuracy effect.

The rejected explicit RHS-direction experiment 1248899 converged in 17
macros / 377 HVPs / 12.85 s, -2181.624464668362 Eh, 481.5 MiB RSS.
It added 29 independent correction directions without reducing total work.
Its method, call site and test were removed, with no retained runtime switch.

## Mathematical prototype

Let the frozen coupled Hessian action be

$$
\mathcal H = \begin{pmatrix} A & B^{\mathsf T} \\ B & C \end{pmatrix}.
$$

These are abstract actions; production does not construct their matrices.
Let $P_\omega$ be the existing signed block-diagonal residual inverse at
the current AH eigenvalue. Store independent orbital directions in $Q$
and their already computed coupling images in $R=BQ$, with
$Q^{\mathsf T}Q=I$. Orthogonalization applies the same linear combination
to the sampled directions and their images; it requires no new $B$ action.
The coupling approximation is

$$
\widehat B=RQ^{\mathsf T},\qquad
U=\begin{pmatrix}Q&0\\0&R\end{pmatrix},\qquad
K=\begin{pmatrix}0&I\\I&0\end{pmatrix}.
$$

It reproduces $B$ on the sampled orbital space. Its transpose is defined
algebraically, ensuring a symmetric approximate coupled operator; it is
not claimed to reproduce every action of $B^{\mathsf T}$.
The signed Woodbury inverse is

$$
\widehat P_\omega
=(P_\omega^{-1}+UKU^{\mathsf T})^{-1}
=P_\omega-P_\omega U
\left(K+U^{\mathsf T}P_\omega U\right)^{-1}
U^{\mathsf T}P_\omega.
$$

This changes only the residual preconditioner. The AH equation, full HVP,
physical metric, nonlinear keyframe gradient, orbital update and accuracy
criteria are unchanged. It introduces no additional coupling/HVP calls and
no history across macro points. Rank is bounded by the existing AH workspace
size and dependence is tested with its existing numerical rank criterion.
Storage is linear in the coordinate counts times that rank, not quadratic
in orbital or structure dimension; no thread-local large replicas are used.

The small signed solve is residual-certified. A singular approximate model
is reported rather than hidden by a different production path. Such a failure
would reject the prototype; no extra empirical spectral regularization has
been introduced. The prototype first needs a dense independent unit reference
for positive and indefinite diagonals, followed by a new FeCl2 short run.
Speedup and use on porphine remain unproven until those measurements succeed.

Prototype job 1248922 completed its bounded 32-macro test without convergence:
650 HVPs, 23.69 s, 495.2 MiB RSS, -2181.624459511211 Eh. Its preconditioner
scope cost 2.558 s. Five unit tests, including the independent sampled-model
inverse reference, passed. This certifies the implemented inverse formula,
not the spectral effectiveness of its approximate model. Compared to the
retained 16/334 result, both convergence and cost deteriorated. The prototype
class, implementation, source registration, driver hooks and unit experiment
were deleted. No failed algorithm or runtime toggle remains in production.
The existing 32-step bound applies only to this short diagnostic input, not
to the production optimizer or the independently running porphine jobs.

The result rejects this particular combination of exact sampled coupling
with an independent diagonal base. It does not refute low-rank response
representations generally or prove a specific spectral pathology. A next
proposal needs a coherent approximation of the diagonal and coupling actions
on the same space, with measured model error, rather than an empirical damping
parameter. Neither the full operator nor its tolerances was changed here.

## Coherent coarse-space alternative to investigate next

The failed construction suggests using **all** components of the already
evaluated actions on a block-separated sample space, not only the coupling.
For every combined direction $(p,q)$, the current HVP separately evaluates
$Ap$, $Bp$, $B^{\mathsf T}q$ and $Cq$. With independent sampled orbital and
structure bases, let $V=\operatorname{diag}(Q_o,Q_s)$ and retain $HV$ by
applying the same orthogonalization combinations to those component images.
These images are available without new Hessian or coupling evaluations.
Existing metric images must likewise be reused, not reevaluated.

For the shifted action $L=H-\omega G$, form only the small coherent operator
$E=V^{\mathsf T}LV$. A symmetric balancing-deflation preconditioner is

$$
P_2=VE^{-1}V^{\mathsf T}
+(I-VE^{-1}V^{\mathsf T}L)P_\omega
 (I-LVE^{-1}V^{\mathsf T}).
$$

Its application needs two small solves and only one base-preconditioner
action. It reproduces the inverse on the sampled coarse space when $E$ is
invertible, without replacing exact coarse diagonal blocks by approximate
diagonals. Unlike using the combined Davidson space alone, the separated
space can contain new directions outside the existing combined space.
At this stage it was an unimplemented proposal, not an improvement claim. First verify
the identity, symmetry, signed shift and rank behavior against independent
dense models, then run a new bounded FeCl2 test. Do not silently regularize
singular $E$ or add damping parameters to force a passing result.

After removing both prototypes, restored-source verification 1248925 passed
all five unit tests and reproduced 16 macros / 334 HVPs, energy
-2181.624464776005 Eh, 11.43 s GNU wall, 487.7 MiB RSS. Its 0.39 s difference
from 1248881 is not an algorithmic speedup: the solving policy is restored
and action counts are identical. The production changes retained at the end
of that round were opt-in scalar diagnostics and automatic result parsing only.

## Implemented coherent coarse inverse

`optimization/preconditioners/coarse.cpp` now implements the balancing formula.
The CIAH image observer is called after an admitted direction's exact HVP;
it receives the already computed metric image. The driver supplies the pure
orbital and pure structure images evaluated within that HVP. No HVP, metric
or coupling replay is added. Every macro constructs a fresh coarse workspace,
bounded by twice the existing AH space capacity; the three shared image
arrays are not duplicated per thread. Rank and symmetry checks use the
existing dependence and declared operator-accuracy contracts, not new
optimization thresholds. The shifted coarse LU solves are residual-certified
and do not discard negative eigenvalues.

Independent dense tests check the balancing matrix, symmetry, exact sampled
inverse identity, signed/shifted operators, dependence rejection and one
base inverse call per application. A CIAH observer test confirms one observer
per admitted HVP and the matching cached metric image. All five existing
compute-node regression tests passed in 1248942 and 1248947.

| Variant | Job | Macros | HVP | Energy / Eh | GNU wall / s | RSS / MiB |
|---|---|---:|---:|---:|---:|---:|
| Retained diagonal reference | 1248925 | 16 | 334 | -2181.624464776005 | 11.43 | 487.7 |
| Coherent coarse prototype | 1248942 | 13 | 257 | -2181.624464862785 | 10.47 | 512.8 |
| Coarse inverse with symmetry certification | 1248947 | 13 | 257 | -2181.624464861511 | 10.44 | 520.6 |

No old executable was rerun. All use 32 bound Hanhai22 CPUs and the original
guess, exact integrals, original CIAH policies, g-infinity tolerance 1e-3 and
energy-change tolerance 1e-7. Nodes differ; wall times are exploratory
single-run comparisons. The last coarse run has 23.1% fewer HVPs and 18.8%
fewer macros; the 8.7% wall improvement is not a controlled scaling study.
The first coarse prototype's preconditioner cost was 0.1906 s. Its recorded
HVP time includes coarse sample ingestion. This is not a hardware-profiler
CPU percentage or an isolated integral-kernel speedup.

### Energy references: preserve both

The existing PySCF default AH endpoint remains -2181.6244647028957 Eh, making
the final coarse endpoint lower by 1.5862e-7 Eh. Do not silently replace it.
The older independent PySCF coupled-Newton population audit already converged
to -2181.624464883733 Eh under the same 1e-3/1e-7 stopping tolerances, using
CAS(10 electrons,8 orbitals), spin=4, Cartesian def2-TZVP, no symmetry, the
same Fe--Cl geometry and the same original input orbitals. Its initial energy
-2181.6176319285287 agrees with the AH reference -2181.617631928641 to
1.123e-10 Eh. The coupled-Newton endpoint is recorded in
`results/diagnostics/fecl2-pyscf-reference-20261011/charges-result/charges.csv`;
the script `charges.py` records macro=24 for a converged final solve (a
nonconverged final solve would be labelled max_cycle_macro=100).

The new final coarse energy differs from that coupled-Newton endpoint by
2.2222e-8 Eh. No rerun, changed threshold or extra convergence criterion was
used to obtain this comparison. The two old PySCF endpoints themselves
differ by 1.8084e-7 Eh; the energy-change criterion is not a certificate of
absolute error relative to the mathematical minimum.

### Other systems and porphine gate

Frozen certified binary 1248947 was used by 1248951's independent array jobs:

| Input | Macros | HVP | Energy / Eh | GNU wall / s | RSS / MiB |
|---|---:|---:|---:|---:|---:|
| F2 | 3 | 23 | -198.751155830526 | 0.32 | 17.4 |
| 241 | 3 | 21 | -230.720590392836 | 3.20 | 841.3 |
| MnF2 | 10 | 137 | -1348.893353223951 | 5.68 | 527.4 |

All converged. These are current-build results, not matched improvements over
an old version; historical runs have different preconditioners and timing
scopes and must be labelled accordingly. A bounded four-macro original-guess
porphine array 1248959 now tests singlet/triplet with the same frozen binary.
No replacement unbounded porphine run is launched before those results are
examined. Existing original-guess long runs remain unchanged.

The full goal remains incomplete: FeCl2 still requires 257 HVPs and neither
porphine NEO energy agreement nor efficient convergence has been verified.
