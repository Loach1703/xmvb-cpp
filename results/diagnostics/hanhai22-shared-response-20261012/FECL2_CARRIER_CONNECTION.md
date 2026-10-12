# Carrier connection in the coupled nonorthogonal model

## Status and purpose

This is a derived diagnostic hypothesis, not an established implementation
bug or an implemented optimizer change. The FeCl2 tail residual is mostly
orbital; neither orbital-block Jacobi nor metric-orthogonal AH vectors
reduced the work. The next test should separate physical integral coupling
from coupling caused by changing the representation of the CI coefficients.
No orbital direction, spin space, retraction or threshold should be changed
before this decomposition is measured and checked.

## Structure and spin-string spaces

At orbital coefficients $b$, let $F(b)$ map VB structure coefficients
$\mathbf C$ into the globally orthogonal spin-string product carrier:

$$
y=F(b)\mathbf C,\qquad M=F^TF,\qquad H=F^T h F.
$$

Here $h$ is the orthogonal-carrier Hamiltonian, $M$ and $H$ are structure
overlap and Hamiltonian matrices, and $y$ is a flattened spin-string product
tensor. These formulas do not introduce a determinant intermediate in the
implementation. Neither $F$ nor a complete dense $H$ or $M$ is proposed as
a production allocation.

For FeCl2, $U_\alpha=\binom87=8$ and $U_\beta=\binom83=56$ give 448
carrier entries. The complete target-spin structure space has dimension
$448-\binom88\binom82=420$, not 448; the extra carrier spin sector must not
be mistaken for an additional variational freedom.

## Exact mixed derivative

Write $r=(h-E)y$ at $\mathbf C^TM\mathbf C=1$. For an orbital direction
$p$, denote $\delta F=DF[p]$ and $\delta h=Dh[p]$. Before projecting onto
the structure tangent, the mixed derivative of the coefficient gradient is

$$
Bp=2\delta F^Tr+
2F^T(h-E)\delta F\mathbf C+
2F^T\delta h\,y-
2\delta E F^Ty.
$$

The last term vanishes after projection onto an overlap-horizontal
structure tangent, because its columns $Z$ satisfy $Z^TM\mathbf C=0$.
In a complete invariant target-spin space, $r=0$ at the relaxed eigenstate.
Then only two terms remain:

$$
Z^TBp=2Z^TF^T(h-E)\delta F\mathbf C+
2Z^TF^T\delta h\,y.
$$

The first is a representation (carrier connection) term, and the second
contains the physical integral response. For an arbitrary selected VB
subspace, $F^Tr=0$ does **not** imply $r=0$. The $\delta F^Tr$ term and any
closure defect must therefore be retained. Completeness must not be inferred
merely from complete unique-spin-string lists.

## Connection direction and its closure certificate

Remove normalization from the basis-induced carrier displacement:

$$
d_y=(I-yy^T)\delta F\mathbf C.
$$

For the unscaled tangent map $F_Z=FZ$, find its metric least-squares
coordinates from

$$
(F_Z^TF_Z)\gamma=F_Z^Td_y,\qquad
e_{\mathrm{closure}}=d_y-F_Z\gamma.
$$

This is an overlap-metric solve, **not** a structure-Hessian response solve.
When the selected structure space is invariant under the spin-scalar
one-body generator and complete in the target spin, the closure defect is
zero in exact arithmetic. Otherwise only a projected connection is obtained;
the remaining carrier direction is physical and must not be discarded.
In code, tangent coordinates are scaled by $\sqrt{2w}$; this factor must be
included explicitly when relating $\gamma$ to `StructureTangentOperator`.

If $K=2F_Z^T(h-E)F_Z$ is the unscaled structure curvature, completeness gives

$$
Z^TBp=K\gamma+2Z^TF^T\delta h\,y.
$$

Thus the coupled forcing can be large even when the physical integral
coupling is modest. This is a possibility to measure, not the conclusion
of the current tests.

## A possible preconditioner, not a new variational space

If the direction-dependent connection defines $\gamma=\Gamma p$, use the
linear shear

$$
\begin{pmatrix}p\\q\end{pmatrix}
=J\begin{pmatrix}p\\s\end{pmatrix},\qquad
J=\begin{pmatrix}I&0\\-\Gamma&I\end{pmatrix}.
$$

The transformed Hessian is a congruence, not the deletion of a coupling:

$$
J^T\begin{pmatrix}A&B^T\\B&K\end{pmatrix}J
=\begin{pmatrix}
A-B^T\Gamma-\Gamma^TB+\Gamma^TK\Gamma & B^T-\Gamma^TK\\
B-K\Gamma & K
\end{pmatrix}.
$$

The full AH metric must transform as $J^TGJ$ as well. Alternatively this
shear may be used only to construct a preconditioner, leaving the real AH
operator and metric unchanged. A block-diagonal metric in the new coordinates
cannot be substituted without acknowledging that it changes the AH model.
No performance benefit is assumed: applying $\Gamma$ and its adjoint may
cost enough to outweigh better conditioning.

## Bounded implementation and measured result

Hanhai22 job 1249390 compiled the diagnostic and passed `test_ciah`,
`test_coupled_structure`, `test_orthogonal_direct_ci`,
`test_sparse_orbital_quotient` and `test_hessian_diagonal` on the compute node.
The new unit checks both a stationary complete-space root and a normalized
nonstationary mixture. It verifies unchanged H/M images, the exact forcing
decomposition, complete-space connection closure and $K\gamma$ at the stationary
root. The nonstationary test explicitly retains the left carrier-residual term.

The implementation reuses the existing generator and sigma results. No extra
orbital HVP or AO-integral traversal is made. Full trace alone captures three
additional carrier/structure blocks with a conservative 16 MiB workspace guard;
normal/scalar trace does not capture or solve this diagnostic. All workers share
the accepted data; there is no per-thread carrier copy. Connection metric solves
use the existing certified square-root-of-machine-epsilon residual rule.
Additional CI and overlap actions in this audit are diagnostic overhead, not
production acceleration. The capped diagnostic is not a new approximation or
response-admission policy.

All 281 admitted orbital directions were captured. The maximum reconstruction
error was $2.1435\times10^{-14}$ in the independent structure coordinates. The
largest certified relative connection-metric residual was
$1.4893\times10^{-8}$. The largest absolute relative closure-squared discrepancy
was $2.1407\times10^{-10}$; small signed negative values are retained in the raw
logs, not clipped to hide metric-solve error.

The following entries are medians of **per-direction ratios** in each macro,
not ratios of averaged vectors. Cross cosine uses the complete inverse metric.

| Macro | Directions | $\|B_{\rm rep}p\|/\|Bp\|$ | $\|B_{\rm int}p\|/\|Bp\|$ | Cross cosine | $\|B_{\rm rep}p-K\gamma\|/\|B_{\rm rep}p\|$ |
|---|---:|---:|---:|---:|---:|
| 11 | 22 | 0.512 | 1.118 | -0.503 | 9.56e-8 |
| 12 | 23 | 0.484 | 1.037 | -0.264 | 1.03e-7 |
| 13 | 26 | 0.549 | 0.945 | -0.309 | 1.25e-7 |
| 14 | 29 | 0.544 | 1.001 | -0.325 | 8.72e-8 |

Basis motion is real and sizable, but does not dominate the total forcing or
produce enormous hidden cancellation. These data do **not** support deleting
it or attributing the slow tail primarily to a spurious coefficient
re-representation. They also do not establish that a full connection shear is
worth its inverse/adjoint costs. No shear or new preconditioner has been enabled.

| Run | Mode | Macros | HVP | Energy / Eh | GNU wall / s | GNU RSS / MiB |
|---|---|---:|---:|---:|---:|---:|
| 1249390, anode18 | Full connection audit | 14 | 281 | -2181.624464869989 | 11.81 | 549.3 |
| 1249395, anode31 | Scalar, final guarded source | 14 | 281 | -2181.624464869471 | 10.55 | 512.3 |
| Archived 1249221 | Scalar, before this diagnostic | 14 | 281 | -2181.624464869146 | 9.93 | 510.0 |

Both new jobs converge with the unchanged 1e-3 gradient / 1e-7 Eh energy
criteria; both pass the five unit executables. The normal-mode energy differs
from the retained PySCF coupled-Newton endpoint by +1.4262e-8 Eh. Timings come
from different nodes and different trace modes; the diagnostic run must not be
used as a performance baseline. The normal run shows no action-count gain and
is not a claimed speedup. No archived optimizer was rerun and no iron-porphine
long run was submitted.

Raw data: `fecl2-connection-1249390/native/{run.err,neo.tsv,output.out,resource.time}`
and `fecl2-ritz-cost-1249395/native/`. Reproduce the direction table with
`summarize_connection.py LOG --output CSV`; the archived `connection.csv` reports
every macro, including inverse certification and action accounting.

The remaining target is the orbital complementary-space residual and how its
expansion/preconditioning represents the physical integral response. The
connection diagnostic narrows that target; it does not solve it.

## Existing code and next bounded diagnostic

`StructureAction::apply_integral_direction` already forms the relative
upper-triangular orbital-transform direction, then
`sigma.apply_one_body_generator(spin_vectors, relative_orbital_direction)`.
This is the basis-induced $\delta F\mathbf C$, distinct from a physical CI
response. Its sigma action and the integral-direction sigma action are
currently added before projection back to structures. The decomposition can
therefore be observed at this point without an additional orbital HVP or
AO-integral traversal.

The diagnostic should report, in the complete structure dual norm:

1. Total coupled forcing, basis-induced forcing and physical integral forcing.
2. Their cross pairing; separate norms alone can hide strong cancellation.
3. The metric solve's unchanged true residual and the carrier closure defect.
4. Reconstruction of total forcing from these terms, with the current
   normalization, state weight and tangent-coordinate factors.

Use FeCl2 first and existing accepted-point data. Do not infer a bug from
large components alone, form a dense coupling matrix, extend the selected
space, or start another iron-porphine long run before this check.
