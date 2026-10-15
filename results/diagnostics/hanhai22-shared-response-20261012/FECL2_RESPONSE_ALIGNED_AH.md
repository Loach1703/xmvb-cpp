# FeCl2: response-aligned AH companion directions

## Motivation and scope

The orbital/structure diagonal and correct local orbital blocks have not
resolved the slow complementary residual. Existing sampled Schur spectra
show up to 99% orbital-curvature cancellation by structure relaxation.
The retained AH trial direction is a combined (p,q) vector. Its q comes
from the previous residual preconditioner and need not approximate the
CI response to the newly admitted orbital p.

This experiment adds that response as a separate AH direction. It is not
the rejected split-Ritz trial, which merely freed the already present q;
not sampled-coupling Woodbury; not a carrier-connection shear; and not a
new diagonal inverse. It does not change the exact coupled action, chart,
physical metric, retraction, keyframe tests or convergence criteria.

## Response-aligned subspace

At one frozen accepted point, denote the coupled actions by

$$
H=\begin{pmatrix}A&B^T\\B&C\end{pmatrix},\qquad
G=\operatorname{diag}(G_o,G_s).
$$

They are operators, never full production matrices. The normal primary
HVP already returns Ap, Bp, B^Tq and Cq separately. Reuse Bp and solve

$$
(C+\mu G_s)s=-Bp,\qquad \mu=\max(0,-\omega),
$$

where omega is the previous AH root, not a new adjustable level shift.
For the usual negative root, C+mu G_s is precisely the shifted structure
block in the AH residual equation. The ground-state horizontal C is
positive semidefinite; strict positivity depends on a gap or positive mu.
The initial zero-shift solve therefore requires a nonsingular horizontal
ground-state block. No fallback for degenerate/indefinite cases is proposed.

The new trial vector and its exact images are

$$
w=\begin{pmatrix}0\\s\end{pmatrix},\qquad
Hw=\begin{pmatrix}B^Ts\\Cs\end{pmatrix},\qquad
Gw=\begin{pmatrix}0\\G_ss\end{pmatrix}.
$$

Together, the existing primary direction (p,q) and w permit an independent
coefficient for the specific induced structure response. A richer trial
space can represent relaxed soft directions more directly; it does not
guarantee better nonlinear convergence. No Schur complement is inserted
into the actual HVP, and no structure coordinate is removed.

The response solve uses the existing matrix-free SPD CG and its certified
true relative residual sqrt(machine epsilon). This is a direction generator,
not an AH stopping criterion. The AH residual is still evaluated using the
complete actual H and G. Approximate response directions would still give
a valid Ritz space when their H/G images are exact; this trial does not
introduce an additional response-accuracy policy or loosen any certificate.

## Exact images without unstable image subtraction

The companion's direct source and images remain together. It is not
orthogonalized by subtracting nearly dependent Hessian images. A separate
Euclidean QR span is retained only to generate independent new primary
directions. The Ritz basis may thus be nonorthogonal, which the existing
generalized AH eigensolver handles with the actual projected metric.
Companion dependence checks retain the existing local and projected-metric
rank rules. Both orientations of H/G projections are still checked;
the symmetry certificate is not widened.

The original count of primary orbital HVP expansions remains at most 30.
Each can supply at most one additional CI direction, so the memory capacity
is bounded by twice that existing count. The first usable-step and microstep
policies still count primary HVP admissions, not the enlarged vector count.
The AH budget is not increased to make a weak primary expansion look better.

## Costs and accounting

- No additional A/B orbital HVP or AO-integral traversal is required.
- Each companion requires a structure CG solve and one B^T action. The
  latter is a structure-response orbital adjoint, not another A action.
- The CG true-residual certification supplies the final Cs/Gs images and
  adjoint inputs. Those are reused, not recomputed by another H/S action.
- The existing primary coarse observer remains unchanged; CI companions
  do not pretend that old component images are their own images.
- Additional primary-span/basis/image arrays are shared, O(d m), with
  m<=60; there is no full C, B, Hessian, integral replication, or per-thread
  large array. The worst-case additional tall storage must still be reported.
- `ciah_response_companion_cost` records CI actions, CG iterations, the
  true residual, shift, and wall scope separately. `orbital_hvp_actions`
  retains actual orbital HVP meaning; more CI work must not be hidden as
  an HVP speedup. The ordinary response counter includes the extra adjoint.

The saved final CG coupling result and a few coordinate vectors are shared
linear-size workspace, not a stored inverse. Structure CG may be very costly
for large active spaces. Thus FeCl2 is the first gate; iron is not run merely
because the linear algebra tests pass.

## Gates and provenance

Independent dense tests use a four-dimensional coupled operator with a
nonidentity metric. They check full H/G images and residuals, the exact
full generalized AH root after two primary actions plus two companions,
unchanged primary action budget, no replay at exhaustion, and rejection of
a dependent companion. Existing five unit executables gate chemistry.

First build/run: Hanhai22 Slurm **1249778**, original FeCl2 guess, exact
integrals, 32 bound cores, same analytic diagonal and outer criteria.
Archived 1249596 (anode31, 12/251, 8.30 s, 516.23 MiB,
-2181.624464871899 Eh) is the comparison; it is not rerun.
Strict PySCF reference remains -2181.624464883733 Eh. Adoption requires
correct endpoint energy and a real wall-time/action improvement, not just
a bigger basis or a solved structure equation.

## Completed result: rejected

Job 1249778 completed on anode31 with all five unit executables passed.
The original FeCl2 guess **converged in 32 macros / 620 orbital HVPs**,
at -2181.624464870668 Eh, with GNU wall **25.10 s** and peak RSS
**543.43 MiB**. Application exit 0 and the convergence message confirm
this was a converged endpoint at the last diagnostic-bound macro, not
an iteration-cap failure. Energy is 1.3065e-8 Eh above the strict reference.
It is much slower than the archived 12/251/8.30 s baseline.

620 companion constructions required **30,744 CI actions**, with median
47 CG iterations and maximum true relative residual 1.48948e-8. Their
explicit construction scope cost **4.2364 s**. The ordinary response
counter is 1,240; it does not include those 30,744 internal CI actions,
which must be reported separately. Small AH spectrum/image reconstruction
cost 0.12316 s. Its increase is not the dominant failed-trial cost.

The macro exit distribution was 23 micro-budget exits, five gradient-growth
exits, and four keyframe-mismatch exits. The last macro still ended at its
budget, with full residual 1.06848e-4. Only 0.082% of that residual squared
was structural. Accurate structure response and a larger exact-image trial
space did not cure the orbital tail. The sampled structure spectrum's
moderate range in older diagnostics is not a full-space condition bound;
47-iteration CG behavior must not be dismissed based on that small sample.

All companion production APIs, response solves, enlarged workspace, altered
admission accounting, and trial tests were removed. No switch or fallback
remains. The physical operators, original AH space and numerical controls
are restored. This test does not justify another structure-response model
or an iron long run.

The next distinct check is an actual **origin coupled HVP versus a centered
difference of the complete coupled gradient**, along the first usable AH
step. The old keyframe audits compared local nonlinear curvature to the
frozen origin model; that difference alone cannot certify the origin HVP.
The new opt-in audit keeps both orbitals and structure coefficients in the
same origin chart, uses the existing two finite-difference spacings, adds
no HVP, and makes no optimizer decision. It is not another fixed-point
energy comparison. Diagnostic evaluation time is not native performance.

## Completed origin derivative audit

Hanhai22 job **1249786**, anode31, used the restored algorithm plus the
opt-in audit. All five unit executables passed. Four accepted macros produced
eight centered-difference comparisons (two spacings per macro). Application
exit 2 is the intended four-macro diagnostic cap, **not convergence**; the
completed batch status is not a convergence certificate. This audited run
is not a timing comparison. Its source, input and executable hashes are
saved with the raw logs in `fecl2-origin-hvp-1249786/native`.

The comparison uses the same fixed orbital lift and structure tangent map
for both sides:

$$
e(\delta)=\frac{g(\delta d)-g(-\delta d)}{2\delta}-H d,
\qquad d^T G d=1.
$$

Here g is the complete coupled gradient pulled back to the immutable origin
chart, not a reoptimized-CI orbital gradient. Report the complete physical
dual norm, so a small coordinate norm cannot conceal a metric-amplified error.

| Macro | Relative HVP error, delta=1e-3 | Relative HVP error, delta=3e-4 | Orbital absolute error, delta=3e-4 | Structure absolute error, delta=3e-4 |
|---|---:|---:|---:|---:|
| 1 | 1.80799e-6 | 1.62854e-7 | 2.75255e-6 | 8.24385e-8 |
| 2 | 3.96734e-7 | 1.91193e-7 | 1.80381e-8 | 8.85899e-10 |
| 3 | 3.09525e-7 | 2.11117e-7 | 2.57242e-8 | 2.88789e-9 |
| 4 | 7.06784e-7 | 1.22001e-7 | 2.59495e-8 | 5.02144e-9 |

Macro 1's full absolute error drops from 3.05722e-5 to 2.75378e-6,
a factor 0.09007, consistent with the squared spacing ratio 0.09.
The other macros reach an absolute error floor of roughly 2e-8. The smaller
spacing's energy-directional-derivative errors span 2.55e-9 to 4.82e-8 Eh.
No systematic origin HVP inconsistency is demonstrated along these actual
AH directions. This is a directional test, **not a proof for every direction**.
It also does not imply that frozen-origin curvature remains exact after a
finite nonlinear keyframe displacement.

Only the opt-in audit and parsing support are retained. It makes no model,
preconditioner, threshold, step or stopping change. The prior parser fields
were checked against HEAD on archived 1249619 and are identical.

## What the micro-budget actually means

The archived retained run 1249619 ends macro 11 at 24 usable microsteps,
26 AH directions, and a full residual 2.86069e-4; the existing micro-gradient
target is 9e-5. The orbital residual is 2.85940e-4 and the structure residual
8.61168e-6: 99.91% of its squared norm is orbital. No microstep is clipped.
The last keyframe correction projects only about 3% of its squared dual
norm into the AH space. Macro 12 then exits on the actual micro-gradient
test, not the budget. Thus "always exits at the budget" is not literally
true; repeated earlier budget exits do show insufficient residual reduction.

PySCF's `mcscf/newton_casscf.py:update_orb_ci` also checks the microstep
budget before its keyframe trigger. Its `soscf/ciah.py:davidson_cc` retains
xs/ax, changes the gradient row of the projected AH problem, and adds only
one new image per unconverged expansion. Our retained implementation follows
these ordering and image-reuse choices. A changed keyframe RHS does not,
by itself, invalidate all frozen-Hessian images. PySCF's Euclidean norm is
appropriate to its orthonormal orbital/CI coordinates; substituting it for
our physical nonorthogonal dual norm would change the problem's metric.
References: [Newton CASSCF source](https://github.com/pyscf/pyscf/blob/master/pyscf/mcscf/newton_casscf.py),
[CIAH Davidson source](https://github.com/pyscf/pyscf/blob/master/pyscf/soscf/ciah.py).

The native 1249619 scopes give orbital dual-norm inversions 0.23973 s,
structure inversions 0.08921 s, AH spectrum/image reconstruction 0.01977 s,
and projection updates 0.02421 s, against 7.40 s total. The structure metric
solve needs at most one CG iteration, with true certified relative residual
8.84e-16. These scopes are not all disjoint from higher-level timers; do not
sum nested scopes into an invented wall total. Neither dense small AH
eigensolving nor full norm inversion explains most of the runtime or the
orbital residual that remains outside the trial space. The remaining target
is effective orbital-space enrichment/model reuse, not weaker certification
or more structure response solves. No iron calculation is authorized by this
negative test as an improvement gate.
