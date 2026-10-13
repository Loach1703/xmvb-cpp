# Origin-model globalization of a nonlinear CIAH sequence

## The rejected-ray defect

The co-iterative augmented Hessian (CIAH) iteration retains curvature images
at one accepted orbital/structure point and refreshes the true gradient at
nonlinear keyframes. Let $\mathbf p$ be the accumulated displacement in that
point's immutable nonredundant chart, and let $\mathbf g_0$ be its gradient.
For the ground-state energy after structure relaxation,

$$
E(t\mathbf p)=E_0+t\mathbf g_0^T\mathbf p+O(t^2).
$$

Here structure coefficients are stationary at the macro origin, up to the
unchanged CI solve accuracy. Orbital/structure gradient and Hessian images
are expressed in the same chart. Finite keyframe gradients need not agree
with the frozen origin prediction.

Descent of an individual microstep relative to its refreshed gradient does
not imply $\mathbf g_0^T\mathbf p<0$ for their sum. If the accumulated trial
is rejected and $\mathbf g_0^T\mathbf p>0$, repeatedly halving that same ray
cannot recover local descent. If its quadratic origin term is also positive,
the origin quadratic model predicts an increase at every positive scale.
This is a globalization defect, not evidence that the origin HVP is wrong.

Saved porphine traces show this failure mode. They must not be confused
with a new failure reproduction or proof of convergence after correction.
PySCF's keyframe procedure does not impose our separate monotonic macro
energy veto; copying its approximate microsteps does not by itself supply
a descending ray for that additional veto.

## Reuse the retained second-order space

Let $Q$ contain the retained Euclidean-orthonormal AH directions, and let
$\mathcal H_0$ and $\mathcal G_0$ denote the matrix-free coupled Hessian and
positive coordinate-metric operators at the origin. Their images are already
available. Define only small projected matrices and a projected gradient:

$$
H_Q=Q^T\mathcal H_0Q,\qquad
G_Q=Q^T\mathcal G_0Q,\qquad
\mathbf g_Q=Q^T\mathbf g_0.
$$

On a failed macro ray, solve the origin augmented problem in this existing
space:

$$
\begin{pmatrix}
0&\mathbf g_Q^T\\
\mathbf g_Q&H_Q
\end{pmatrix}
\begin{pmatrix}\alpha\\\mathbf z\end{pmatrix}
=\lambda
\begin{pmatrix}1&0\\0&G_Q\end{pmatrix}
\begin{pmatrix}\alpha\\\mathbf z\end{pmatrix}.
$$

Select its lowest root. For the raw tail $\mathbf u=Q\mathbf z$, the first
row gives

$$
\mathbf g_0^T\mathbf u=\lambda\alpha.
$$

When $\mathbf g_Q\ne0$, the augmented matrix has a negative Rayleigh
quotient, hence its lowest root is negative. Orient the tail against the
origin gradient. A nonzero reference component then gives first-order
descent. A root with $\alpha=0$ instead satisfies

$$
\mathbf g_0^T\mathbf u=0,\qquad
\mathbf u^T\mathcal H_0\mathbf u
=\lambda\mathbf z^TG_Q\mathbf z<0,
$$

and is a genuine negative-curvature direction. It need not capture the full
gradient or solve the Newton equation.

The usual microstep reference-weight filter cannot certify descent: it may
discard a weak-reference negative root and select a positive root. In one
dimension, $H_Q=-100$, $G_Q=1$, and $g_Q=0.001$ give a negative root with
$|\alpha|\approx10^{-5}$ and a positive root with $|\alpha|\approx1$.
Consequently the restart uses the lowest root without that filter; the
ordinary microstep root-selection policy is unchanged.

## Stable clipping and true-energy acceptance

Let $L(\mathbf u)$ be the existing maximum orbital/structure component
measure, and let $d_{\max}$ be its unchanged microstep cap. Use

$$
\tau=\min\left(\frac{1}{|\alpha|},
                    \frac{d_{\max}}{L(\mathbf u)}\right),\qquad
\mathbf p=\tau\mathbf u,
$$

with $\tau=d_{\max}/L(\mathbf u)$ when $\alpha=0$. Compare
$L(\mathbf u)$ with $d_{\max}|\alpha|$ before dividing by $\alpha$.
This produces the same component-clipped AH step for a nonzero reference
weight without constructing an enormous intermediate vector.

The image uses the same cached combination,
$\mathcal H_0\mathbf p=\tau(\mathcal H_0Q)\mathbf z$.
Restore the origin-model prediction, the orbital/structure displacement,
and the recycled last direction together. Do not retain the previous
keyframe energy anchor for this different ray.

The normal microstep candidate and the lowest evaluated descending
keyframe are tried first. The cached origin restart is used at most once
when those fail and the accumulated ray is non-descending, the keyframe energy
model has failed as described below, or when its
ordinary backtracking has been exhausted. Zero accumulated steps can also
use this ray. It introduces no additional HVP, metric action, physical
metric inverse, preconditioner call, or accuracy threshold. It requires one
small augmented eigensolve and linear combinations of existing images.

The restarted ray is still tested with a fully relaxed CI energy under the
original true-energy acceptance rule. No model prediction is treated as an
accepted energy, and no extra outer convergence condition is introduced.
The recorded origin Ritz residual uses full cached images, including the
component outside $Q$, but belongs to the unscaled augmented eigenvector.
It is not a stationarity certificate for the clipped accepted macro step;
the latter is explicitly recorded as approximate.

If the retained space contains neither a nonzero gradient projection nor
negative curvature, it cannot provide an origin ray without enrichment.
The method does not silently claim convergence in that case. Likewise,
accurate trial CI relaxation remains necessary; this globalization repair
does not establish lowest-root ordering from a CI residual alone.

## Independent checks

Unit tests compare cached images/residuals with dense small-matrix references,
cover weak-reference and zero-reference negative roots, and assert that
operator, preconditioner and physical-norm callback counts do not increase.
For the nonlinear example

$$
f(x)=x-\frac{10}{3}x^2-900x^3,
$$

a cumulative ray $p=0.03$ is uphill at every positive scale up to one,
whereas the clipped descending origin direction $p=-0.03$ lowers energy.
This is a ray/globalization witness, not a claim that the molecular
keyframe scheduler produces a prescribed artificial microstep sequence.

## Use the already evaluated keyframe energy

Let $\mathbf p_a$ be the most recent keyframe displacement, with true energy
$E_a$ and true pulled-back gradient $\mathbf g_a$. The retained curvature
predicts the next keyframe energy at $\mathbf p$ as

$$
\widehat E(\mathbf p)=E_a+\mathbf g_a^T(\mathbf p-\mathbf p_a)
+\tfrac12(\mathbf p-\mathbf p_a)^T\mathcal H_0(\mathbf p-\mathbf p_a).
$$

The previous gradient-agreement test alone could permit another gradient
refresh when $\widehat E(\mathbf p)<E_a$ but the already computed true
keyframe energy satisfies $E(\mathbf p)\ge E_a$. In that case stop the
approximate microstep sequence before refreshing its energy anchor. Retain
the current displacement and its image for the ordinary relaxed-CI trial.
The original gradient-mismatch rollback takes precedence when it applies.

This sign test has no empirical parameter and requires no new physical
evaluation. It is a conservative globalization decision, not a proof that
further microsteps could never recover. In particular, keyframe energies use
the current coupled structure coefficients, whereas the subsequent macro
trial relaxes them. A rising keyframe energy therefore does not imply that
the relaxed trial must be rejected. If it is rejected, try the evaluated
descending checkpoint and then the cached origin ray before backtracking.

In the saved triplet late-point diagnostic, a keyframe predicted a decrease
of $1.119$ mEh but its true energy rose by $3.778$ mEh. Gradient agreement
still passed the existing test. Five subsequent approximate microsteps
extrapolated an $83.7$ mEh macro decrease, and the ensuing trial CI solve
remained busy when the diagnostic was cancelled. The sign safeguard stops
at the informative keyframe; it is not an extra outer convergence condition.

## CI residual preconditioning in the physical carrier

Let $F$ map structure coefficients into the orthogonal spin-string-product
carrier. The structure overlap is $M=F^TF$. For a complete selected spin
sector, the existing occupation-block representation applies $M^{-1}$
without forming a global dense overlap. Let $D_E$ be the existing
regularized carrier Hamiltonian diagonal shifted by the CI Ritz energy
$E$. The CI residual is a structure covector,

$$
\mathbf r=(H-EM)\mathbf C.
$$

Its preconditioned correction is

$$
\delta\mathbf C=-P_E\mathbf r,\qquad
P_E=M^{-1}F^T D_E^{-1}F M^{-1}.
$$

Both Riesz maps are needed. Using only $F^TD_E^{-1}F$ as an inverse would
confuse structure vectors with covectors. The carrier diagonal need not
preserve a pure-spin subspace, so $P_E$ is a preconditioner, not an exact
inverse of $H-EM$. Projection through $F^T$ returns a valid structure-space
correction. All Davidson actions, residual tolerances, root selection and
physical residual checks remain unchanged.

This representation is admitted only for complete spin sectors on the
direct-CI path, with the existing factor and workspace byte limits.
Arbitrary selected subspaces retain their general structure-coordinate
preconditioner. Occupation blocks and exterior transformations are shared
by CI and the accepted-point model. The transient action has stable owned
storage; its metric operator is rebound once when the action is moved into
the accepted-point context. No factor arrays or integrals are copied.

Tests compare $P_E$ against an independently assembled small-matrix
reference, check action moves and borrowed lifetimes, and check generic
Davidson callback signs and failure handling. Molecular regression also
covers state averaging and TNHVP, because they move the forward context
through different ownership paths.

## Computational validation (14 October 2026)

Hanhai22, 32 OpenMP threads, BLAS single-threaded, exact integrals,
projected-gradient tolerance $10^{-3}$ and energy-change tolerance
$10^{-7}$ Eh. CI solve tolerances and the 256 MiB exact-exchange memory
budget are unchanged. Timings below are single runs, not repeated samples.

The three-macro porphine diagnostic starts from each spin's saved late-point
Molden orbitals, **not** the original input guess. Job 1252645 used frozen
build 1252634, before the forward-action ownership correction. That
ownership correction was subsequently required by the F2-SA/TNHVP paths;
it does not change the numerical kernels used in this bounded test.

| Spin | Initial energy / Eh | Energy after three steps / Eh | HVP | Wall / s | Peak RSS / GiB |
|---|---:|---:|---:|---:|---:|
| Singlet | -2244.731433756452 | -2244.733688467149 | 35 | 192.37 | 38.62 |
| Triplet | -2244.780021402979 | -2244.782339737575 | 44 | 210.50 | 38.63 |

Both tests ended at the deliberately imposed diagnostic `max_iterations=3`;
neither is a converged calculation. The singlet control on the same node,
with matching input/basis hashes, took 287.28 s for three steps/35 HVP and
ended at -2244.733688440202 Eh. The 26.9 nEh endpoint difference is smaller
than the unchanged energy-change tolerance. The triplet control did not
finish its first trial before cancellation after 37 min 19 s; its total
speedup cannot be reported as a completed paired timing.

The triplet initial CI solve decreased from 308.8747 to 0.6392 s. The new
metric factors occupy 3,067,680 bytes in 14,025 occupation blocks, with a
maximum block dimension of 28. This is a phase-level result, not a 483-fold
speedup of the entire calculation. Five fully relaxed triplet macro trials
in the new three-step test took 13.93 s in total. Its energy-sign safeguard
activated once; the cached origin-ray restart was not needed. Thus these
molecular data validate the energy-sign stop and CI preconditioning, not
activation of every restart branch. Restart branches have independent
small-matrix unit coverage.

Final frozen build 1252661 passed ten unit executables and the PySCF F2
reference CTest on a Slurm compute node. All eight native regressions in
job 1252666 converged, including F2 state averaging, TNHVP and block-LBFGS.
The six paired NEO energies differ from the saved current-code reference
by less than $10^{-9}$ Eh. Performance is not uniformly improved: FeCl2
changed from 12/251 steps/HVP to 13/276, while MnF2 changed from 13/189 to
13/160. The complete paired records are in
`benchmarks/orbital-preconditioner/ciah-ci-regression-20261014.csv`; the
bounded porphine records are in `ciah-late-point-20261014.csv`.

Job 1252674 runs both spins from the original input guesses with frozen
build 1252661. It has no three- or thirty-macro diagnostic cap (the existing
application safety limit is unchanged). Final convergence, energy agreement
with the saved same-guess PySCF results, and actual origin-restart activation
remain to be established from these runs. The superseded energy-sign-only
triplet diagnostic 1252612_1 was cancelled after 26 min 04 s: it had
accepted two steps but still spent minutes in subsequent unpreconditioned
CI trials. No unrelated Slurm jobs were cancelled.
