# FeCl2: incremental AH projections and the remaining micro-budget exits

## Scope and reference

Hanhai22 job 1249512, anode31, 32 bound cores. Original FeCl2 input,
exact integrals, analytic Hessian diagonal, unchanged nonorthogonal chart,
Hessian action, keyframe policies and convergence thresholds. All five unit
executables passed. Only the new version was run; 1249478 is archived data
from the same node. No iron-porphine calculation was submitted.

The installed PySCF 2.14.0 `soscf/ciah.py:davidson_cc` and
`mcscf/newton_casscf.py:update_orb_ci` were inspected. In particular:

- The micro-budget test precedes the keyframe test in PySCF too.
- PySCF retains its trial vectors and Hessian images after RHS changes.
- PySCF updates only the new row/column of its projected operator.

Source reference: [PySCF CIAH](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/soscf/ciah.py)
and [coupled Newton CASSCF](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/mcscf/newton_casscf.py).
This comparison does not import orthogonal orbital rotations into VBSCF.

## Equivalent projection update

For the frozen accepted-point operators, store basis columns and images

$$
Q_m=(q_1,\ldots,q_m),\qquad Y_m=H Q_m,\qquad Z_m=G Q_m.
$$

The original code rebuilt $Q_m^T Y_m$ and $Q_m^T Z_m$ and copied all three
tall matrices at every Ritz solve. Admission of $q_{m+1}$ changes only the
new row and column. The new implementation evaluates both orientations,
for example

$$
(H_Q)_{i,m+1}=q_i^T Hq_{m+1},\qquad
(H_Q)_{m+1,i}=q_{m+1}^T Hq_i.
$$

The two orientations are not assumed equal: the original symmetry certificate
still checks them. Positive dot sums used in its roundoff bound are cached
alongside the projection entries. No certificate tolerance is changed.
The current RHS still updates $Q_m^T g$ on every solve.
Step and image reconstruction use the stored vector columns directly.

For full coupled dimension $d$ and $m$ admitted columns, cumulative tall
projection work changes from $O(d m^3)$ to $O(d m^2)$. Each new admission
needs $O(dm)$ projection work. The additional four bounded small matrices
occupy 28,800 bytes at the existing 30-column cap, shared by the workspace;
there are no thread-local tall-matrix replicas. This changes neither the
AH equation nor the number/accuracy of required operator images.

The new unit test independently rebuilds the complete projected generalized
AH problem after each admission and RHS refresh, compares its selected Ritz
step and full-space residual, verifies image reuse, and rejects an explicitly
asymmetric operator. Existing nonidentity-metric tests remain active.

## Native results

| Version | Macros | HVPs | Final energy / Eh | GNU wall / s | GNU process peak RSS / MiB |
|---|---:|---:|---:|---:|---:|
| Archived 1249478 | 14 | 282 | -2181.624464869630 | 9.93 | 492.59 |
| Incremental 1249512 | 14 | 280 | -2181.624464867196 | 9.92 | 507.45 |

The archived HVP count is **282**, independently summed from `neo.tsv` and
the admitted-direction events. The earlier narrative recorded 281; it was
a reporting mistake, not a benchmark rerun or a different convergence rule.
The energy difference is $2.434\times10^{-9}$ Eh. Different dot-product and
vector-accumulation ordering does not produce a bitwise-identical trajectory;
the independent AH test, actual endpoint and unchanged accuracy contract are
the correctness evidence. No claim of total runtime or RSS improvement is made.

| Accumulated stage / s | Archived | Incremental |
|---|---:|---:|
| Projection update + Ritz projection/certificate | 0.16159 | 0.04081 |
| Small AH spectrum and step/image reconstruction | 0.01778 | 0.02341 |
| Complete AH residual/model norm | 0.25934 | 0.25910 |
| Structure dual norm, included in norm consumers | 0.11340 | 0.11310 |
| Admitted-direction HVP + coarse observer | 2.48473 | 2.50891 |

The new projection-update timer is included in the first line; it is not
hidden outside the old Ritz timer. Timers are nested: do not add structure
norm again to complete norm cost. Process RSS above excludes batch compiler
RSS; Slurm reported zero MaxRSS for the application step and is not used as
the process-memory measurement.

## Why the budget still binds

This update deliberately leaves the micro-budget unchanged. Macros 12–14
still end there. Archived macro 14 ended with AH residual
$1.2370\times10^{-4}$, above the existing micro target $9\times10^{-5}$;
the orbital contribution was $1.2091\times10^{-4}$ and the structure
contribution $2.6144\times10^{-5}$. The limiting residual is mostly orbital,
not an expensive incomplete structure metric inverse.

The small AH Galerkin solve itself is accurate: the archived maximal
relative projected residual is $3.26\times10^{-12}$. Its small eigensolver
and structure norm together are inexpensive. Thus budget exits cannot be
explained by spending the time inside either of these two solvers.

The sampled Schur audit already shows strong orbital–structure cancellation
and a broad orbital/relaxed spectral spread. The keyframe audit shows most
true-gradient corrections lie outside the current AH space. These establish
where the iteration struggles, but do not prove a defective gradient or a
particular preconditioner fix. Fresh sample-image reconstruction and local
rank-two secant updates have been tested and rejected; they must not be
reintroduced as a supposed new solution.

The next unresolved mathematical question is alignment of the analytic
fixed-orthogonal-carrier orbital model with the raw VB coefficient coupled
coordinates. A representation-connection congruence needs its forward and
adjoint maps, normalization/scaling and actual projected metric derived and
independently validated before an optimizer change. Merely increasing the
micro cap, relaxing a residual tolerance or replacing the norm is not the
repair supported by these measurements.
