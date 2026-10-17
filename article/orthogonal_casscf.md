# Orthogonal CASSCF: shared CIAH and direct-CI model

## Scope and ownership

The initial implementation is real, exact-integral, single-state, complete
active space (CAS). It reuses `mcscf/ci`, the public CIAH Ritz workspace and
microstep driver, the common Davidson solver, AO/active integral responses,
and the analytic fixed-state rotation diagonal. It does not reuse a VB
objective or construct VB structures. State averaging is outside this
initial implementation. The shared `.xmi` frontend selects this model with
the `CASSCF` flag instead of `VBSCF`; core/active spans of the supplied
orbital seed are preserved while completing an AO-metric orthonormal frame.

`casscf/orbitals` defines the rotation chart. `casscf/derivatives` supplies
the physical energy, gradient and coupled actions. `casscf/driver` performs
the model's CI relaxation and physical macro acceptance. The immutable AO
operator is borrowed; active integral/cache storage is shared by actions at
one orbital point, never copied per OpenMP worker. The accepted-point exact
2e cache borrows owned active coefficients, so its owner is nonmovable.

## Orbital and CI coordinates

Let $U$ be a complete MO frame satisfying $U^TSU=I$, where $S$ is the AO
overlap. Only inactive–active, inactive–virtual and active–virtual rotations
are retained. For their real parameters $p$, define $K(p)^T=-K(p)$ and

$$
U(p)=U\exp K(p).
$$

Internal rotations within each sector are redundant for this complete CAS.
This assumption is not extended to an arbitrary selected VB subspace.

The normalized CI coefficient array $c$ has alpha-string rows and
beta-string columns. Its norm and inner products are Frobenius products.
Let $P_S$ be the exact projector onto the target total spin, not merely a
fixed $M_S$ sector. The origin CI tangent projector is

$$
Qz=P_Sz-c\langle c,P_Sz\rangle.
$$

No dense matrix for $Q$ is formed. For a coefficient direction $q$, the
finite state and its normalization factor are

$$
c(q)=\frac{c+Qq/\sqrt{2}}{N(q)},\qquad
N(q)=\sqrt{1+\tfrac12\|Qq\|^2}.
$$

The $1/\sqrt{2}$ scale makes the CI–CI Hessian action $Q(H_a-E_a)Q$ rather
than twice that operator. Here $H_a$ is the active Hamiltonian, $E_a$ its
expectation, and the total energy also includes inactive and nuclear terms.

## Exact AO responses

With inactive coefficients $U_i$, active coefficients $A$, and the shared
linear symmetric AO response $L(D)=2J(D)-K(D)$,

$$
D_i=U_iU_i^T,\quad F=h+L(D_i),\quad
h_a=A^TFA,\quad
E_i=E_{\mathrm{nuc}}+\operatorname{tr}[D_i(2h+L(D_i))].
$$

The exact two-electron active transform supplies $H_a$. Sigma and its
integral adjoint use the common spin-string connection graph. If
$\gamma=\partial E_a/\partial h_a$ is its symmetric one-electron adjoint
and $\Gamma$ its packed two-electron adjoint, the raw coefficient gradients
are

$$
G_i=4FU_i+2L(A\gamma A^T)U_i,\qquad
G_a=2FA\gamma+G_2(\Gamma).
$$

$G_2$ is the existing exact active-integral pullback. Its directional
derivative and the active-integral variation share the accepted 2e cache.
The variations of $D_i$ and $A\gamma A^T$ are evaluated in one AO-response
batch; the changing-CI adjoint is supplied by the fused sigma/adjoint sweep.
No relaxed orbital response is added on top of the explicit CI tangent.

## Joint Hessian action

Write $g_o$ for the rotation gradient and
$r_c=Q(H_a-E_a)c$. The joint gradient is

$$
g=\begin{pmatrix}g_o\\\sqrt{2}r_c\end{pmatrix}.
$$

At the origin, for $q=Qq$, the coupled action is

$$
H\begin{pmatrix}p\\q\end{pmatrix}
=\begin{pmatrix}
A_op+\sqrt{2}B^Tq\\
\sqrt{2}Q\,\delta H_a[p]c+Q(H_a-E_a)q
\end{pmatrix}.
$$

$A_o$ is the fixed-CI orbital Hessian action, including the exponential
coordinate connection. $B^Tq$ denotes the orbital pullback of
$\langle q,H_ac\rangle$. These are actions, not stored Hessian matrices.
The response sigma and the integral adjoint for $B^Tq$ are fused. The
orbital and CI diagonal preconditioners contain no relaxed CI response.

## Finite keyframe pullback

All microsteps of a macroiteration use its fixed origin coordinates.
For a raw orbital gradient $G$ at $U(p)$ and $B=U^TG$, the rotation
cotangent is obtained from the adjoint Frechet derivative

$$
\nabla_K E=\mathcal L_{\exp}(K^T,B),\qquad
\frac{\partial E}{\partial p_{ij}}
=(\nabla_K E)_{ji}-(\nabla_K E)_{ij}.
$$

One block exponential supplies this entire cotangent. There is no
exponential per rotation coordinate. At the origin its directional
derivative is the skew packing of

$$
U^T\delta G+\tfrac12(K^TB+BK^T).
$$

For the same finite displacement, the CI keyframe gradient in the origin
coordinates is

$$
g_q(p,q)=\frac{\sqrt{2}}{N(q)}Q
\left[H_a(p)c(q)-E_a(p,q)c(q)\right].
$$

Thus both exact keyframe components are derivatives of the same finite
energy, not current-frame gradients inserted into an old chart. The common
CIAH driver controls subspace expansion, clipping, keyframes and microstep
progress. CI is relaxed by the shared Davidson solver for physical macro
acceptance. No extra terminal stability condition is introduced.

## Coupled coarse inverse and initial directions

The orthogonal model uses the same `mcscf::CoarsePreconditioner` as VBSCF.
For an admitted joint direction $x=(p,q)$, one physical action already supplies
the two separate images

$$
H\begin{pmatrix}p\\0\end{pmatrix}
=\begin{pmatrix}A_op\\\sqrt{2}Q\delta H_a[p]c\end{pmatrix},\qquad
H\begin{pmatrix}0\\q\end{pmatrix}
=\begin{pmatrix}\sqrt{2}B^Tq\\Q(H_a-E_a)q\end{pmatrix}.
$$

These split samples enlarge the coarse space beyond
the mixed Ritz directions, with no extra AO response or sigma sweep. Let
$V$ contain its independent samples, $L=H-\lambda I$, and
$Z=V(V^TLV)^{-1}V^T$. For the original diagonal inverse $D^{-1}$, the common
balancing correction applies

$$
P^{-1}=Z+(I-ZL)D^{-1}(I-LZ).
$$

It is symmetric and satisfies $P^{-1}LV=V$ on the retained coarse space.
The basis, images and correction belong to one frozen macro point and are
discarded at the next point. There is no extra physical action, dense Hessian
or per-thread copy of the coarse model.

The retained samples must also support a numerically symmetric projected
operator. QR transforms a source-image error $\Delta Y$ into
$\Delta YR^{-1}$; independence alone does not control that amplification.
Candidate coarse samples therefore satisfy the existing projected symmetry
certificate before admission. The dot-product term uses the full coordinate
dimension through $\gamma_n=nu/(1-nu)$, where $u$ is unit roundoff. Samples
that do not satisfy this certificate are not retained in the optional
preconditioner. The complete physical CIAH action and its separate symmetry
test are unchanged. Certified small Hessian and metric projections are
cached on admission and reused for every subsequent shifted coarse solve.

Only the last accepted microstep is recycled as the next macro's initial
direction, following the existing VB and PySCF `dr0` policy. MO indices are
not canonicalized or reordered. The CI eigenvector sign is aligned with the
finite predictor, and the recycled CI tail is projected into the new point's
target-spin tangent. This is an initial search direction, not transported
curvature, a reused physical gradient or a Hessian image. Every new model
computes its own exact images. Backtracking or a checkpoint restart updates
the recycled step consistently with the actually accepted path.

## Verification contract

Small nonsymmetric test Hamiltonians check energy derivatives, finite
keyframe pullbacks, joint HVP finite differences, Hessian symmetry and
orbital–CI reciprocity. Singlet, unbalanced doublet and maximum-spin cases
check target-spin preservation and removal of the CI scale direction.
Production gradients and HVPs are analytic; finite differences appear only
in these tests. Native Slurm/PySCF comparisons use one binary AO/MO fixture
for both methods and retain their initial/final energies, action counts,
timings, memory and source/environment hashes.
