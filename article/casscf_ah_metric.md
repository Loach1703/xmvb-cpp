# CASSCF augmented-Hessian metric and step scale

## Scope

**Status:** tested as an isolated experiment on 17 October 2026, not adopted
as the production default. The derivation is valid, but the measured
convergence/performance benefit was not consistent across the tested input
paths. The exact tested implementation and extra tests are preserved in
`benchmarks/casscf-metric-20261017/metric.patch`; production code remains
the preceding identity-AH-metric version. Results:
`benchmarks/casscf-metric-20261017/RESULTS.md`.

This note derives the orthogonal complete-CAS limit of the existing VBSCF
inactive/active span metric. It is a metric for orbital subspaces, **not**
the norm of a many-electron wave-function derivative. The experiment changes
the generalized augmented-Hessian (AH) metric only. Energy, gradient, Hessian
actions, finite orbital/CI updates, component clipping, accuracy thresholds,
keyframe policy and outer convergence criteria retain their existing forms.

## Orthogonal orbital frame

Let $S$ be the atomic-orbital overlap matrix and
$U=(U_c,U_a,U_v)$ a complete molecular-orbital frame, with inactive, active
and virtual columns satisfying $U^TSU=I$. An independent real rotation has
the antisymmetric generator

$$
K(p)=\begin{pmatrix}
0&-P_{ac}^T&-P_{vc}^T\\
P_{ac}&0&-P_{va}^T\\
P_{vc}&P_{va}&0
\end{pmatrix}.
$$

Here $P_{ac}$, $P_{vc}$ and $P_{va}$ contain active–inactive,
virtual–inactive and virtual–active rotation amplitudes, respectively;
$p$ packs their independent elements in `RotationChart::pairs()` order.
The finite update remains $U(p)=U\exp K(p)$. At the origin,

$$
\dot U_c=U_aP_{ac}+U_vP_{vc},\qquad
\dot U_a=-U_cP_{ac}^T+U_vP_{va}.
$$

Only the three inter-sector rotations are independent for this complete
CAS. No selected-VB-subspace closure assumption is used.

## Pullback of the span metric

Use any AO whitening factor $W$ with $W^TW=S$, solely to define Euclidean
projectors $\Pi_c=WU_cU_c^TW^T$ and $\Pi_a=WU_aU_a^TW^T$. The existing
VB span metric is the sum of the two Grassmann metrics,

$$
\|p\|_{G_o}^2=\tfrac12\|\dot\Pi_c\|_F^2+
                 \tfrac12\|\dot\Pi_a\|_F^2.
$$

Equivalently, its two horizontal features at this orthogonal frame are
$W\dot U_c$ and $W\dot U_a$, because both variations are perpendicular to
their own subspace. Substitution gives

$$
\tfrac12\|\dot\Pi_c\|_F^2
=\|P_{ac}\|_F^2+\|P_{vc}\|_F^2,
\qquad
\tfrac12\|\dot\Pi_a\|_F^2
=\|P_{ac}\|_F^2+\|P_{va}\|_F^2.
$$

Therefore

$$
\boxed{\|p\|_{G_o}^2=2\|P_{ac}\|_F^2+
\|P_{vc}\|_F^2+\|P_{va}\|_F^2.}
$$

Inactive–active rotations move both subspaces and occur twice; the other
two classes move only one retained subspace. Thus $G_o$ is a diagonal
operator with weights $2,1,1$ in the existing packed rotation coordinates.
There is no numerical whitening, AO projector construction, dense metric
storage or iterative inverse needed to apply this result.

## CI scale and the joint AH equation

Let $c$ denote the normalized, target-spin-pure CI array and $Q_c$ its
target-spin tangent projector. The current CI chart is

$$
c(q)=\frac{c+Q_cq/\sqrt{2}}{\sqrt{1+\|Q_cq\|^2/2}}.
$$

The common coupled convention scales a coefficient displacement by
$\sqrt{2}$, so $\sqrt{2}\,\dot c=Q_cq$ at the origin. The CI metric is
therefore the identity on that tangent; it does not acquire another factor
of two. The joint metric is $G=\operatorname{diag}(G_o,I_c)$. The physical
gradient $g$ and Hessian action $H$ are unchanged, while the AH root obeys

$$
\begin{pmatrix}0&g^T\\g&H\end{pmatrix}
\begin{pmatrix}\alpha\\v\end{pmatrix}
=\omega\begin{pmatrix}1&0\\0&G\end{pmatrix}
\begin{pmatrix}\alpha\\v\end{pmatrix},\qquad
d=v/\alpha.
$$

Equivalently, $(H-\omega G)d=-g$. Every metric image supplied to the
existing CIAH workspace and coarse preconditioner must use this same $G$.
The analytic diagonal preconditioner uses
$h_{ii}-(\omega-\epsilon_{\mathrm{shift}})G_{ii}$, with its existing
$\epsilon_{\mathrm{shift}}=10^{-8}$ and denominator protection unchanged.
This is the matching shifted operator, not a new empirical preconditioner.
The first unshifted seed is unchanged.

## Physical component cap is distinct from metric length

At an orthogonal frame, the VB complete-frame component cap reduces to
the largest actual rotation amplitude, not the largest metric-whitened
coordinate. In the common CI scale, the existing joint cap is

$$
m(d)=\max\left(\max_{ij}|p_{ij}|,
                    \max_{\alpha\beta}|(Q_cq)_{\alpha\beta}|/\sqrt{2}\right).
$$

It remains clipped at the existing $0.03$. In particular, an
inactive–active angle of $0.03$ is not additionally divided by $\sqrt{2}$.
The metric length used by the existing keyframe interval is instead

$$
\ell(d)=\sqrt{d^TGd}.
$$

These measure different properties: the cap controls the largest physical
coordinate increment; the length is the origin metric's quadratic measure
of accumulated subspace movement, not an exact finite geodesic distance.
Conflating them would change two experimental factors.

## Isolated validation

Production changes are limited to packed metric weights, metric images,
the matching diagonal shift and the existing keyframe length callback.
The coordinate residual norm and reported outer gradient norm remain
Euclidean as in the baseline; no physical-metric inverse is introduced.
Gradient/HVP and finite-keyframe derivative tests must remain unchanged.
Independent tests compare the packed metric to explicit projector
derivatives, confirm the CI scale, and show that clipping is unchanged.
FeCl2 comparisons use the stored baseline results, never a rerun of the old
binary. A reduction in iterations or wall time is an empirical outcome,
not a consequence guaranteed by this derivation.
