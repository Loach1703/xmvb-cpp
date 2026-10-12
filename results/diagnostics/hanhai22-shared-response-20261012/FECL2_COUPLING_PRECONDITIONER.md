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
This is an unimplemented proposal, not an improvement claim. First verify
the identity, symmetry, signed shift and rank behavior against independent
dense models, then run a new bounded FeCl2 test. Do not silently regularize
singular $E$ or add damping parameters to force a passing result.

After removing both prototypes, restored-source verification 1248925 passed
all five unit tests and reproduced 16 macros / 334 HVPs, energy
-2181.624464776005 Eh, 11.43 s GNU wall, 487.7 MiB RSS. Its 0.39 s difference
from 1248881 is not an algorithmic speedup: the solving policy is restored
and action counts are identical. The production changes retained in this
round are opt-in scalar diagnostics and automatic result parsing only.
