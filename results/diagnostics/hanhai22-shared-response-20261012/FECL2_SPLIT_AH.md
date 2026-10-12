# FeCl2: why a larger AH trial space did not fix the tail

## Controlled question

The coupled action already computes the two images for an orbital-only
direction and a structure-only direction. Their sum is the image used by
CIAH. Experiment 1249010 used those two images as separate Ritz directions,
without additional integral traversals, HVPs, metric actions, or looser
numerical criteria. All images stayed in the immutable nonorthogonal chart.
The HVP budget and the first usable-step test continued to count evaluations,
not the doubled number of basis vectors.

Independent tests compared the generalized AH root, both operator images,
and the full physical residual with a dense four-dimensional reference.
They also checked a zero sector and exhaustion after two evaluations rather
than after two vectors. These and the four existing regression executables
passed on a Hanhai22 compute node. Passing these tests establishes the linear
algebra, not an improvement in nonlinear optimization.

## Native results

All runs use the original FeCl2 guess, exact def2-TZVP integrals, 32 bound
cores, Davidson, the analytic diagonal, gradient tolerance 1e-3 and energy
tolerance 1e-7 Eh. The two references below are archived; neither was rerun.
Times are exploratory single-run measurements on different nodes.

| Job / implementation | Macros | HVP | Total energy / Eh | GNU wall / s | GNU RSS / MiB |
|---|---:|---:|---:|---:|---:|
| 1248925, diagonal before coherent coarse inverse | 16 | 334 | -2181.624464776005 | 11.43 | 487.7 |
| 1248947, retained coherent coarse inverse | 13 | 257 | -2181.624464861511 | 10.44 | 520.6 |
| 1249010, split Ritz space / base diagonal | 16 | 295 | -2181.624464843057 | 11.07 | 532.1 |

The split-space result differs from the retained result by 1.8454e-8 Eh.
It reduces HVPs relative to the old diagonal-only run, but is inferior to the
retained coarse inverse in macros and HVPs. It is not an accepted optimization.
The experimental splitter, evaluation counters and tests were removed from
production and the execution mirror; no fallback switch remains.

The coarse inverse was intentionally omitted from the split-space trial:
its coarse space would then be the actual Ritz space. On a Galerkin residual,
its inverse correction inside that space is eliminated by the subsequent
orthogonalization. It would add work without enlarging the next complementary
direction in exact arithmetic. This is why the comparison is between two ways
of using the same sampled sector information, not just a bigger basis with
an unrelated arbitrary preconditioner.

## Cost attribution and tail diagnosis

These are explicit wall-clock scopes, not CPU-profiler sample percentages.
AH residual norms are included in the all-structure-norm row; do not add the
overlapping rows to estimate total time.

| Scope / s | Retained 1248947 | Split 1249010 |
|---|---:|---:|
| AH projection and checks | 0.1357 | 0.3113 |
| Spectrum and Ritz image combination | 0.0167 | 0.0389 |
| All structure dual norms | 1.0027 | 1.0421 |
| Residual preconditioning | 0.1912 | 0.1011 |
| Coupled expansion HVP scope | 2.291 | 2.186 |

The retained HVP scope includes coarse-image ingestion; the split scope
includes the image-split callback but not its subsequent image orthogonalization.
These scopes therefore cannot certify an isolated kernel speedup.

The full structure-norm solve is not free, but is not the explanation for
hundreds of HVPs. Projected AH solves are still less costly. Increasing Ritz
rank improved neither the overall trajectory nor the micro-budget exits.

In the retained run, macros 9--11 stop at `keyframe_mismatch`; macros 12 and
13 stop at `micro_budget`. Thus not every tail exit has the same cause.
In the last macro the physical source gradient norm is 5.9818e-4, and the
final predicted norm is 2.4663e-4, still above the unchanged micro stop 9e-5.
There is no component clipping in these final two macros. Raising the budget
or relaxing this stop would not diagnose the failure to contract the residual.

For retained macros 9--13, the orbital contribution accounts for 99.07--99.66%
of the squared physical norm of the keyframe gradient correction. This
decomposition uses the same block-diagonal physical metric as the complete
norm; it is not a comparison of differently scaled Euclidean sectors.
The split trial shows the same orbital-dominated correction in its final
four macros. This is evidence to investigate the frozen orbital curvature
versus the nonlinear keyframe, not proof of an incorrect analytic HVP.

## PySCF and next verification

The installed `pyscf/soscf/ciah.py::davidson_cc` solves the projected AH problem
and retains its operator images after yielding a usable approximate step.
`pyscf/mcscf/newton_casscf.py::update_orb_ci` has both a gradient-dependent
microstep budget and real-gradient keyframes. A budget exit is a macro model
refresh, not a certificate of a solved full Newton equation. Those policies
were not changed in this experiment. PySCF's Euclidean norms in orthonormal
coordinates do not justify dropping the complete physical norm in our chart.

The next bounded diagnostic uses the retained frozen binary and its existing
directional audit at nonzero keyframes. Compare the analytic energy derivative
with centered differences, and the origin HVP with centered differences of
the *current chart gradient*, along the actual resumed direction at two
displacements. This separates derivative inconsistency from genuine curvature
variation before changing model-refresh logic. It is not another fixed-point
energy comparison or an old-version performance rerun. No new porphine job
is warranted by the failed split-space experiment.

## Evidence

Job 1249010, node anode18, exit 0; five unit executables passed.
Binary SHA256: be5d2754f5dc5bc1b6610f2f25f1b5bc5b2a43344a4814c111dd91fbe4521e0f.
Input SHA256: cd89d6f43c4ed3bd1dc221433426636c0723e4a4db8ba42297db88708583233c.
Raw logs, flags, input and native GNU time are retained in
`fecl2-ritz-cost-1249010/`. Scalar accounting is reproducible with
`summarize_ritz_cost.py fecl2-ritz-cost-1249010/native`.

## Nonzero-keyframe diagnostic result

The follow-up 1249025 used the frozen retained binary, not the split experiment.
It finished with 13 macros / 257 HVPs and energy -2181.624464863389 Eh.
Its runtime is diagnostic overhead and is not a performance comparison.
The initial submission 1249020 mistakenly combined scalar tracing with the
expensive-audit flag: that flag requires full tracing and no directional
audit was emitted. It is not derivative evidence. The corrected batch enables
full tracing and fails if the requested directional events are absent.

At the first accepted nonzero keyframe in each late macro:

| Macro | Origin/current HVP difference divided by origin HVP norm | Energy directional derivative error at delta=1e-3 / Eh |
|---|---:|---:|
| 9 | 2.926% | 8.54e-9 |
| 10 | 8.979% | 1.06e-8 |
| 11 | 3.888% | 2.57e-9 |
| 12 | 3.204% | 1.62e-8 |
| 13 | 0.577% | 4.70e-9 |

The HVP differences agree at delta=1e-3 and 3e-4. This is measured curvature
variation between the origin and the current nonlinear point, not evidence
that an origin Hessian must equal the local Hessian. The energy derivatives
support consistency of the current gradient along these actual directions.
This check does not certify all directions or the analytic local HVP itself.

The bounded independently assembled full structure metric and its certified
matrix-free inverse gave the same keyframe decisions. Diagonal versus full
metric diagnostics also produced no keep/trigger disagreements in this run.
Thus neither an inaccurate inverse nor changing the norm certificate explains
these particular budget exits. No production accuracy was changed.

These observations narrow the problem but do not solve it. A frozen model
does change noticeably at intermediate tail keyframes; in the last macro
the observed curvature change is much smaller, yet the complementary residual
still contracts too slowly. Both the orbital preconditioner's spectral
effectiveness and curvature reuse across keyframes need to be distinguished;
the data do not justify declaring a gradient bug or simply increasing budgets.
