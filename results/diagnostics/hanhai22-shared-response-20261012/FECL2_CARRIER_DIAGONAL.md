# Orthogonal-carrier Jacobi inverse in horizontal VB coordinates

## Motivation and scope

The raw VB structure diagonal and the orthogonal-carrier Hamiltonian diagonal
define different approximate inverses. PySCF 2.14.0 builds the latter through
`make_hdiag` in `mcscf/newton_casscf.py`. Its factors of two/state weights
cannot simply be copied: our structure coordinates already contain
$\sqrt{2w_s}$. The present change replaces only the eligible structure
preconditioner. The exact coupled Hessian, physical metric, nonorthogonal
orbital chart, keyframe/energy acceptance and all tolerances are unchanged.
The rejected carrier-connection shear is not included.

Reference: [PySCF coupled Newton CASSCF](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/mcscf/newton_casscf.py).
The normalized CI diagonal also has a residual-dependent term in PySCF;
this work uses an approximate inverse at a relaxed accepted CI point, not
a claim that a carrier Jacobi model is the exact coupled Hessian diagonal.

## Symmetric pullback

At a fixed accepted orbital point, let $F$ map VB structure coefficients
directly to the orthogonal unique-spin-product carrier. Let $Z$ be the
existing independent horizontal structure basis and define

$$
L=FZ,\qquad G_s=L^TL.
$$

Here $G_s$ is the physical horizontal structure metric, not its diagonal
approximation. For a structure coordinate covector $r$, the minimal-norm
carrier representative satisfying $L^T\rho=r$ is

$$
\rho=L G_s^{-1}r.
$$

For the carrier diagonal $D=\operatorname{diag}(H_{\mathrm{orth}})-E_s-\omega$,
where $E_s$ is the active-space state eigenvalue and $\omega$ is the current
AH eigenvalue shift, pull the carrier inverse back as

$$
P_s r=G_s^{-1}L^T D^{-1}L G_s^{-1}r.
$$

Both Riesz maps are necessary. This expression is symmetric for a signed
diagonal too. It reduces to $(L^TDL)^{-1}$ when $L$ is square and invertible.
For a proper carrier subspace it is an approximate inverse, not the inverse
of the projected diagonal model. Its output remains inside the prescribed
horizontal VB space; no spin/structure direction is added. In particular,
an intermediate carrier diagonal need not commute with the spin operator:
the exact adjoint/Riesz map projects its output back into the allowed sector.

Each $G_s^{-1}$ application uses the existing factorized preconditioner and
the same complete true-residual-certified metric solve. No precision is
relaxed. Eligibility requires an admitted complete spin sector and a bounded
carrier workspace; arbitrary selected spaces retain their existing raw
structure model, rather than pretending their metric is the identity.

## Diagonal without sigma probes

For alpha/beta occupied orbital sets $a,b$, the scalar carrier diagonal is

$$
d_{ab}=\sum_{p\in a}h_{pp}+\sum_{p\in b}h_{pp}
+\frac12\sum_{\sigma\in\{a,b\}}\sum_{p,q\in\sigma}
\big[(pp|qq)-(pq|qp)\big]
+\sum_{p\in a,q\in b}(pp|qq).
$$

With string occupation matrices $N_\alpha,N_\beta$, Coulomb matrix
$J_{pq}=(pp|qq)$, and same-spin energies $e_\alpha,e_\beta$, compute

$$
d=e_\alpha\mathbf1^T+\mathbf1e_\beta^T+N_\alpha JN_\beta^T.
$$

This uses existing occupation lists and accepted orthogonal active integrals.
The shared scalar carrier diagonal is built once per accepted action, not
once per preconditioner call. There is no complete H/S matrix, determinant-pair
matrix, extra orbital HVP or electron-integral transformation. The direct
spin-string-to-structure contraction remains unchanged. Only one carrier
column is processed at a time. Admission counts the shared diagonal, carrier
workspaces and occupation/GEMM temporaries against the existing action budget,
including unequal alpha/beta string counts. There are no thread-local replicas.
The existing signed denominator protection is retained unchanged.

## Independent validation

Hanhai22 1249543 and build-only 1249552: all five unit executables passed.

- Occupation-only diagonals are compared elementwise with independent
  Hamiltonian sigma probes in the small direct-CI fixture.
- The nonorthogonal two-orbital fixture independently constructs its full
  carrier transformation and Hamiltonian reference, then checks the entire
  pulled-back preconditioner matrix and symmetry.
- Both single-state and equal-weight two-state horizontal constraints are
  tested. Removing either Riesz map would fail these comparisons.

The final memory-admission accounting changes only eligibility, not the admitted
model or arithmetic. It was built and unit-tested in 1249552 before submitting
the regression and iron-porphine arrays with that frozen executable.

## FeCl2 native benchmark

Both runs: anode31, 32 bound cores, same original input and exact integrals.
1249512 is archived data, not a repeated old-version run.

| Structure inverse | Macros | HVPs | Final energy / Eh | GNU wall / s | GNU process peak RSS / MiB |
|---|---:|---:|---:|---:|---:|
| Raw structure diagonal, 1249512 | 14 | 280 | -2181.624464867196 | 9.92 | 507.45 |
| Carrier Jacobi pullback, 1249543 | 11 | 222 | -2181.624464858105 | 8.76 | 485.12 |

The macro count falls by 21.4%, HVP count by 20.7%, and native wall time by
11.7%. Endpoint difference is $9.091\times10^{-9}$ Eh. The new energy is
$2.563\times10^{-8}$ Eh above the PySCF coupled-Newton reference
(-2181.624464883733 Eh), inside the unchanged energy accuracy contract.
Physical final gradient norms differ (1.58e-4 versus 3.11e-4); both are outcomes
of the same existing projected-gradient and energy-change stopping criteria,
not a new criterion.

Structure inverse application costs increase from 0.217 s to 0.292 s summed
over AH preconditioner calls. The gain comes from fewer complete orbital
actions and accepted-point rebuilds, not a cheaper individual inverse.
All eleven macros still exit at the micro budget: the underlying tail issue
is improved, not cured, and no claim of rapid quadratic convergence is made.

## Other converged regression cases

New array 1249553 uses the frozen 1249552 build; comparison is archived
1249412. Nodes differ, so timings are exploratory, not matched speedup proofs.

| Case | Archived macros/HVP | New macros/HVP | Archived/new wall / s | Archived/new peak RSS / MiB | Absolute energy difference / Eh |
|---|---:|---:|---:|---:|---:|
| F2 | 3/23 | 3/23 | 0.33/0.35 | 17.52/17.77 | <1e-12 |
| 241 | 3/21 | 3/21 | 2.76/2.56 | 829.11/804.91 | 7e-12 |
| MnF2 | 10/144 | 10/137 | 5.00/5.33 | 496.62/521.96 | 2.09e-10 |

All application exits were zero and endpoint energies were unchanged within
numerical precision. MnF2 has fewer actions but is not demonstrated faster.
There is no universal speedup claim and no molecule-specific selector.

## Iron-porphine follow-up

1249554_0 (singlet) and 1249554_1 (triplet) run concurrently on Hanhai22
from the original input orbitals, each bounded to four macros for diagnosis.
The last verified scheduler state is running. These are not long convergence
runs and must not be counted as converged results. Compare with the archived
four-step traces first; keep the original-initial-guess PySCF endpoint/time
references rather than reporting restart step counts as total convergence.
