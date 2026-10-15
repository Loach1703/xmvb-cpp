# FeCl2: exact sampled Hessian with a base-model complement

## Why this is distinct

The retained inverse-balancing model reproduces the exact inverse on all
sampled columns LV, but that property does not determine its extension into
the unknown complement. Before-learning inverse audits showed substantial
error on newly admitted directions. Completing structure response explicitly
did not fix the predominantly orbital residual. This diagnostic changes only
the complementary-space inverse model, not the true coupled action, physical
metric, chart, keyframes, stepsize policy, AH space or stopping tolerances.

It is neither the rejected coupling-only Woodbury inverse nor the rejected
multivariate inverse SR1 model. The former used exact B samples with
uncorrected diagonal A/C. The latter used an inverse-defect least-change
completion. This model completes the **Hessian** on a Euclidean orthogonal
complement; all sampled A, B and C columns remain exact.

## Completion and inverse

Write L=H-omega G for the current-root shifted actual operator, V for the
orthonormal separated orbital/structure sample basis, Y=LV, E=V^T Y,
and Q=I-VV^T. The base inverse is P0, and D=P0^-1 is only an abstract
model; it is never explicitly formed. Define

$$
\widehat L=VEV^T+QYV^T+VY^TQ+QDQ.
$$

For a symmetric sampled action, this has the exact interpolation property

$$
\widehat L V=Y=LV.
$$

Only the unknown complement block is replaced by the original analytic
diagonal model. The off-diagonal samples QY are retained in the full physical
coordinate dimension, rather than discarding their components outside V.
This is the symmetric least-change completion of the base Hessian in the
unknown block, not a claim that the complement is already accurately modeled.

Assuming the constrained blocks are invertible, the complement inverse is

$$
P_Q=P_0-P_0V(V^TP_0V)^{-1}V^TP_0.
$$

It satisfies V^T P_Q=0 and P_Q V=0. Define the small Schur operator

$$
F=E-Y^TP_QY.
$$

Block elimination then gives, for any residual r,

$$
a=F^{-1}(V^Tr-Y^TP_Qr),\qquad
\widehat L^{-1}r=P_Qr+(V-P_QY)a.
$$

The result is symmetric and satisfies the same exact sampled inverse relation
as balancing, but can have a different action outside V. Negative eigenvalues
are retained. Singular constrained Gram or Schur models reject the prototype;
there is no fallback, damped eigenvalue, new rank policy, or empirical shift.
Small solves retain the existing roundoff residual certificate.

## Cost and gate

The prototype reuses V, HV and GV from the frozen accepted point; no HVP,
coupling or integral image is added. At each changed root, however, it needs
2m+1 base inverse calls instead of one, where m<=60 is the already bounded
coarse rank. Structure carrier pullback inverses may make those extra calls
expensive. This cost must be included, not described as free Woodbury reuse.
The two additional tall arrays have shape d-by-m, shared serial workspace;
no full Hessian, metric, integral or per-thread large replica is allocated.

The independent dense unit builds the complete matrix above and directly
inverts it, checking the callback result, symmetry, interpolation, signed
base models, nonidentity physical metrics, changed AH shifts and exact base
call counts. Five existing unit executables gate the original-guess FeCl2
diagnostic on Hanhai22, 32 bound cores. Archived 1249596 (12/251, 8.30 s,
516.23 MiB, -2181.624464871899 Eh) is not rerun. Strict PySCF remains
-2181.624464883733 Eh. This is a testable model, not a speedup claim.
Do not launch iron from unit-test success or merely fewer orbital actions.

## Completed experiment: rejected

Hanhai22 job **1249812** ran on anode31, the same node as archived 1249596.
All five unit executables passed, including the independent signed/shifted
dense completion reference. The original-guess chemistry trial then hit the
existing 32-macro diagnostic bound without convergence: **302 HVPs**,
**-2181.622006640423 Eh**, **16.97 s** GNU wall, **535.04 MiB** GNU RSS.
Final physical gradient norm was 0.0434982. Application exit 2 and the
explicit max-iteration termination are not convergence. The energy is
2.45824 mEh above the strict PySCF endpoint, but this is an unfinished
trajectory, not a different converged energy.

The 32 macros ended with 22 gradient-growth exits and 10 micro-budget exits.
Preconditioning cost 3.29269 s; expansion HVP scope cost 2.58744 s. Complete
orbital/structure norms cost 0.36735/0.11750 s and AH spectrum/image work
0.01523 s. These scopes can nest and are not summed into a wall total.
The additional small-inverse work did not buy effective enrichment.

Exact interpolation and symmetry alone do not ensure a useful unknown-block
completion. In particular, even with positive E and a positive complement
model, E-Y^T P_Q Y can be indefinite. Balancing has a different advantage:
positive E and a positive base inverse give a positive balancing inverse.
This algebra explains a possible failure mechanism, **not a measured proof**
that artificial complement curvature caused every observed growth exit.
No damping, eigenvalue cut, alternative threshold or fallback was tried.

The completion API, inverse, extra tall workspace and altered unit tests
were removed as a whole; the original balancing implementation and tests
are restored. Only this evidence document remains. No iron run was launched.
