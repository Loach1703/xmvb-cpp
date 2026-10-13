# Remove work without consumers from fixed-state keyframes

Base: master 4aedc3c5, retaining the coordinate-norm CIAH policy of 6bcc5bb4.
No AH thresholds, negative-curvature handling, HVP, nonorthogonal chart,
keyframe triggers, or external convergence conditions were changed.

## Necessity before speed

An internal keyframe evaluates a prescribed, generally nonstationary structure
vector. Its consumers need the scalar energy, orbital gradient, normalized
coefficients, and structure residual. They do not diagonalize the structure
Hamiltonian or consume a Davidson preconditioner or eigensystem report.

The change removes that unused work:

- `StructureActionMode::ActionOnly` constructs the same physical operator
  without a diagonal model. The normal Davidson constructor behavior remains.
- Fixed-state gradient preparation does not build a Davidson Jacobi diagonal,
  synthetic root columns, SCF eigenvector copies, or overlap diagonal reports.
- Only state energies are contracted, not an unused full state-by-state matrix.
- Coefficient normalization scales the already computed Hamiltonian and overlap
  images. A second action on normalized columns is mathematically unnecessary.

There is one action per fixed-state keyframe and no new per-thread large array.
The accepted-point eigensolver, output, and derivative preparation are unchanged.
This is a removal of unused work, not an approximate-gradient algorithm.

## Is the gradient refresh itself necessary?

Saved production-policy scalar trace: `ciah-coordinate-audit-1250554/QOCGUM-audit`.
The first macro's second refresh has predicted gradient norm 0.182804 and actual
norm 0.748142; mismatch 0.727186 stops further microsteps. The third macro uses
eight internal refreshes; at its end the frozen-origin model gradient norm is
0.040375 while the refreshed prediction is 0.000263. The last macro's refresh
changes the gradient by only 1.567e-7.

These observations distinguish substantial nonlinear updates from low-impact
updates. They do not prove how a trajectory without a refresh would converge,
nor provide an a priori rule for skipping an expensive refresh. Consequently,
no keyframe is removed and no new empirical threshold is introduced here.
An outer relaxed energy/gradient evaluation and a fixed-C internal keyframe
are different operations and cannot simply be substituted for one another.

The `neo.tsv` keyframes counter includes the final accepted-point gradient once
per macro. Internal fixed-state calls therefore equal that count minus macros
in these native runs. QOCGUM's count 17 means 13 internal calls plus 4 accepted
outer gradients, not 17 internal evaluations.

## Validation

Hanhai22 build/unit job 1250646; six executables passed:
`test_keyframe`, `test_ciah`, `test_coupled_structure`,
`test_orthogonal_direct_ci`, `test_sparse_orbital_quotient`,
`test_hessian_diagonal`.

New tests cover action-only complete and incomplete spin carriers against the
existing dense/reference actions; column-normalization equivalence; fixed-state
energy/residual scaling; nonstationary orbital and structure finite differences
for two weighted states with nonconsecutive labels. No CI stationarity assumption
is used in these derivative tests. Synthetic tests are not new molecular
fixed-point comparisons and are not a claim about a new state-average optimizer.

Native array 1250648, indices 1,4,9,11,12. All five inputs and basis hashes match
the saved reference array 1250538. Exact Libcint integrals, 32 OpenMP threads,
single-thread BLAS, NEO, Davidson, analytic Hessian-diagonal preconditioner,
gradient tolerance 1e-3 and energy tolerance 1e-7 Eh. Old binaries were not rerun.
All runs converged; timings are end-to-end application wall times. RSS is GNU
time process maximum, in MiB. Nodes use AMD EPYC 9654 but are not matched nodes;
single-run timing differences include runtime/node variation.

| Case | Steps old/new | HVP old/new | Wall s old/new | RSS MiB old/new | New energy Eh |
|---|---:|---:|---:|---:|---:|
| FeCl2 | 12/12 | 254/252 | 8.963/7.261 | 549.9/518.7 | -2181.624464871145 |
| 241 | 3/3 | 21/21 | 2.272/2.154 | 792.4/808.5 | -230.720590392829 |
| QOCGUM | 4/4 | 52/52 | 15.230/12.789 | 1922.5/1957.5 | -480.094766681859 |
| CERRAS | 5/5 | 51/51 | 34.043/34.135 | 2026.7/1941.2 | -397.093435004110 |
| LOFLEA | 4/4 | 34/34 | 81.751/74.030 | 2696.2/2705.9 | -399.106646001531 |

Maximum endpoint difference is 3.51e-9 Eh. FeCl2 changes by two HVPs and one
internal keyframe; do not describe that timing change as a pure kernel speedup.
No universal RSS reduction is established. Historical Hanhai25 timings are in
the CSV as a separate comparison, not a matched timing baseline; its historical
FeCl2 endpoint is different and not a correct-energy speed target.

## Remaining keyframe cost

Separate scalar-only new-binary audit 1250649_15 (QOCGUM): 13 internal keyframes,
13 combined H/S actions, zero diagonal storage. Total keyframe time is 1.1472 s,
including 0.97575 s active-space work. Within that active scope, operator setup
is 0.04271 s, H/S action 0.03784 s, and active adjoint 0.13624 s. Nested timers
must not be summed. These timers do not include all macro-level forward work.

Keyframe runtime is roughly 9% of this audit's end-to-end runtime. Even removing
all keyframe work at unchanged trajectory could only save that portion; actual
removal changes the trajectory, so this is a cost bound, not a predicted gain.
There is no evidence to justify another round of small keyframe-action tuning.
Any subsequent change to refresh frequency should address nonlinear-model reuse
and total action counts, and first demonstrate that its removed refreshes are
not needed. The safe cleanup in this report is retained independently of that
future algorithm decision.

Data: `keyframe-lean-comparison.csv`, `keyframe-lean-audit.json`,
`keyframe-saved-audit.json`, `keyframe-lean-accounting.txt`. Parsing is automated
by `compare_lean.py` and `summarize_keyframes.py`; raw outputs remain in the
job-labelled directories and the remote frozen build directory.
