# FeCl2: verified carrier connection, rejected inverse congruence

## Question

The analytic orbital diagonal is evaluated with orthogonal active-carrier
density matrices, whereas the coupled VB equation uses raw structure
coefficients. Could aligning these representations in the preconditioner
remove the excessive AH direction extensions?

This is a preconditioner hypothesis, not a demonstrated defect in the exact
Hessian, gradient or orbital update. It was tested with the original input,
exact integrals and unchanged CIAH policies/accuracy on Hanhai22.

## Representation map

Let $F(b)$ map raw VB structure coefficients to the orthogonal unique-spin
product carrier, and $Y=F(b)C$. At the accepted point, let $Z$ be the existing
horizontal structure coordinate basis and $G_s=Z^T F^T F Z$. For the
single equally weighted state, the scaled structure coordinates use
$\delta C=Zq/\sqrt{2}$. An orbital tangent $p$ changes the carrier through
$\delta F(p)C$, even when raw $C$ is fixed. Its horizontal representation is

$$
\Gamma p=G_s^{-1}Z^TF^T\sqrt{2}\,\delta F(p)C.
$$

The factor $\sqrt{2}$ is the existing coupled-coordinate normalization,
not an adjustable parameter. The trial required an admitted complete spin
sector, so the horizontal carrier change can be represented exactly. Merely
complete alpha/beta string lists would not prove structure-space closure.

For active overlap $O=R^TR$, the upper-triangular generator is

$$
K=\delta R R^{-1},\qquad
K_{ij}=(R^{-T}\delta O R^{-1})_{ij}\ (i<j),\qquad
K_{ii}=\tfrac12(R^{-T}\delta O R^{-1})_{ii}.
$$

The carrier derivative is the existing one-body spin-string generator
$\delta Y=E(K)Y$. There are at most $n_{\mathrm{act}}(n_{\mathrm{act}}+1)/2$
channels. Forward generator contractions and the accepted orbital-preparation
adjoint supply factors $W,V$ with $\Gamma=WV^T$. No electronic-integral,
Hamiltonian sigma, complete coupling or orbital HVP evaluation is needed to
build them. The same certified physical Riesz inverse is used for $G_s$.
The FeCl2 factors had 36 channels and occupied 1,157,184 bytes, shared by
the accepted-point workspace. Allocation admission included the temporary
spin-product carrier; no thread-local copies or full coupled matrices existed.

## Symmetric trial inverse

To hold the horizontal carrier fixed rather than raw coefficients, use

$$
\begin{pmatrix}p\\q_{\rm raw}\end{pmatrix}
=T\begin{pmatrix}p\\q_{\rm carrier}\end{pmatrix},\qquad
T=\begin{pmatrix}I&0\\-\Gamma&I\end{pmatrix}.
$$

Only the approximate inverse was replaced by the congruence

$$
P_{\rm raw}=T P_{\rm base} T^T.
$$

This preserves symmetry (including a signed AH inverse). It is *not* an
exact Hessian inverse: integral-dependent coupling remains, and the full
physical metric under the shear is not block diagonal. The actual AH equation,
metric, HVP, keyframe gradient and orbital retraction were left unchanged.
No secant fitting, Schur-response approximation or threshold change occurred.

## Validation

- All five unit executables passed in 1249527.
- An independent explicit $5\times5$ signed inverse checks $TPT^T$ and its
  symmetry without relying on the factorized implementation.
- Finite differences of a two-orbital Cholesky transform check the energy-free
  generator overlap against the existing full integral-direction connection,
  for stationary and nonstationary coefficients.
- Four-macro full-trace job 1249532 compares the factored $\Gamma p$ against
  the connection extracted from each actual full coupled action. All 55
  comparisons passed numerically; maximal physical relative difference was
  $7.574\times10^{-15}$. This diagnoses the sampled directions, not every
  possible tangent. The four-step run stopped at its diagnostic iteration cap
  (application exit 2); it is not a converged benchmark.

## Native convergence result

Both runs used anode31, 32 bound cores, identical original orbitals, analytic
orbital diagonal and the same energy/gradient criteria. The original result
is archived 1249512, not a repeated old-version calculation.

| Inverse | Macros | HVPs | Final energy / Eh | GNU wall / s | GNU process peak RSS / MiB |
|---|---:|---:|---:|---:|---:|
| Retained base, 1249512 | 14 | 280 | -2181.624464867196 | 9.92 | 507.45 |
| Carrier congruence, 1249527 | 22 | 280 | -2181.624464861861 | 12.57 | 510.33 |

The energy difference is $5.335\times10^{-9}$ Eh. The same total number of
orbital actions is spread over more macros, so accepted-point rebuilds cost
more. Congruence-factor setup itself totalled 0.411 s; preconditioner
applications were 0.221 s versus 0.217 s. Analytic-diagonal builds increased
from 1.936 s to 2.928 s.

The stop distribution also changed:

| Micro stop | Retained base | Carrier congruence |
|---|---:|---:|
| Gradient growth | 1 | 13 |
| Keyframe mismatch | 2 | 0 |
| Micro budget | 11 | 9 |

The last three macros of the trial still hit the micro budget. Reducing its
frequency elsewhere is not an improvement when it is replaced by early
gradient-growth exits and extra outer iterations.

## Decision and remaining question

The connection reconstruction is correct on the tested directions, but this
particular approximate inverse is worse. All trial APIs, the congruence class,
driver changes and their now-unused tests were deleted. The local and cluster
source were restored to c001ee8d's retained algorithm, and the cluster executable
was restored from the frozen 1249512 artifact. No baseline workload was rerun.
Unrelated user RI edits were preserved. No iron-porphine long run was submitted.

This rules out *this representation alignment alone* as a practical cure. It
does not prove that the analytic diagonal is exact in raw coordinates, nor
that all possible connection-aware preconditioners must fail.

A separate unresolved difference from PySCF is the structure preconditioner:
PySCF uses the Hamiltonian diagonal of the orthogonal spin-product carrier;
we use raw VB structure diagonal expectations. These are different approximate
inverses, even though their exact coupled equations represent the same physical
space. Before changing it, derive the adjoint/Riesz maps needed to pull a carrier
diagonal inverse into horizontal structure coordinates; do not copy a carrier
diagonal elementwise into the raw VB coordinate basis.
