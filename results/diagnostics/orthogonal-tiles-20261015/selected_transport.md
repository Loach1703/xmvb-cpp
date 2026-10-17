# Selected-space CI transport without sector completion

Parent: `2e83399a`. This change removes an optional preconditioner, not the
selected-space projection, CI transport, orbital derivatives, or stopping
tolerances. All builds and numerical tests run through Hanhai22 Slurm.
Historical executables are not rerun. Native timing/RSS and scalar phase
accounting serve different purposes and are recorded separately.

## Necessary algebra

Let F contain only the selected VB columns in the common orthogonal
representation used by finite CI transport. It is generally rectangular;
the selected span need not be invariant under exterior transformations.
The selected Gram is M = F^T F. For the transport target z, solve

    M c = b,       b = F^T z.

Differentiating that same equation and its reverse gives

    M dc = db - (dM)c,
    M^T lambda = c_bar,
    c_bar^T dc = lambda^T db - lambda^T (dM)c.

These are selected-space equations. They require neither unselected VB
columns nor a complete spin-product numerical array. The existing overlap
actions and adjoints still evaluate the right-hand sides. They do not imply
that every action costs only O(n_selected), or that transformed support is
small for every selection.

Completing the columns to a spin-pure basis gives a Gram with blocks

    M_full = [ M   K; K^T   L ].

Its principal inverse is (M - K L^-1 K^T)^-1, not M^-1 in general. The deleted
RestrictedMetricInverse built this completed-sector inverse action and used
its selected principal block as a positive preconditioner. It therefore did
not eliminate the original selected solve. Its complement generation,
full-product metadata, zero-extended right-hand sides, exterior inverse
passes and full-product projection are not necessary transport operations.

## Implementation boundary

`CiTransport::solve` retains the exact complete-selected-space factorization
when already available. Otherwise it applies PCG to the selected M action
using the positive accepted-point overlap diagonal. That diagonal remains an
SPD preconditioner at a trial chart; it need not equal the trial diagonal.
The existing true-residual stopping target is unchanged. There is no
completed-sector preconditioner, budget admission or alternate transport.

PCG coefficient vectors have selected dimension. Its cost is the selected
M action cost times the measured iteration count; removing preconditioner
setup does not guarantee lower total runtime. In particular, nearly full
selected spaces may require more iterations. This must be measured, not
hidden by relaxing the solve or by silently accepting a projected residual.

The finite transport already provides a known coefficient representative
c_base = c_origin + dc_independent. Start the original solve from c_base,
not zero. This is algebraically the correction equation

    M delta_c = b - M c_base,       c = c_base + delta_c.

The existing solver certifies the full original residual and retains its
original RHS-relative stopping target; it does not redefine the tolerance
relative to the smaller correction RHS. At zero independent CI increment,
the target and F_new*c_base agree through first orbital order. With a CI
increment, the mixed orbital/CI defect is retained too. This is an initial
solution, not an assumption that a selected subspace closes or that a defect
can be discarded. No numerical parallelism/amplitude test is introduced.

This change does not remove the complete-space carrier Hessian-diagonal
preconditioner, the exact complete-space metric inverse, or broad-support
accepted/orbital-direction arrays. They are separate remaining work.

## Validation contract

`test_selected_metric` replaces tests of the deleted completion. It checks
selected M actions and solves against an independent small dense Gram solve,
including one structure, singlets, unequal alpha/beta occupations, rescaled
columns, occupation-mixed columns and a changed trial chart with the old
positive diagonal. A separate small full-space oracle retains the proof
that a principal full inverse is not the restricted inverse.
The existing finite-transport tests check a nonclosed selected span,
normal-equation projection, energy-gradient differences, HVP differences
and multiple roots. Spin-impure columns are used only to test the general
Gram solver, not to admit a different physical spin sector in production.

Native 240 COV/ION/FULL and the complete-space FeCl2 control are compared
with saved final-layout runs1261850/1261856. FeCl2 controls the unchanged
complete-space path; it does not exercise the deleted restricted completion.
Use `tests/mcscf/selected_transport.sbatch` for build/tests and independent
native run jobs, with dependency on the successful build.

## Initial deletion-only experiment

Job1261985 passes all11 tests. Native jobs1261986/1261987 remove completion
but still start selected solves from zero. The five240 cases converge with
unchanged steps/HVP and energy differences below3e-12 Eh versus saved1261850.
Their SCF times are7.012,9.345,11.803,12.802,9.297 s for COV,ION01,ION02,ION03,
FULL. The four truncated spaces are5--10 percent slower than saved6.685,
8.660,10.769,11.733 s; FULL is lower than saved9.741 s although its mathematical
path is unchanged. These are single native measurements on the same shared
anode26, not proof of per-function attribution. Peaks remain1516--1532 MiB;
the full molecular-process memory peak is not materially reduced.

Separate scalar-phase job1261992 records6148 and7046 selected M actions and
1.318/1.668 s transport solve wall for ION02/ION03. Older saved phase capture
1260422 records940/688 actions and0.625/0.494 s with completed preconditioning.
That older capture is not the immediate parent, but directly shows the
preconditioner's iteration-reduction role; deletion alone is not a speedup.
Neither old executable is rerun.

## Final seeded implementation

Job1262005 passes all11 focused tests, including independent minor-based Gram
actions and finite selected-space transport gradient/HVP/multistate checks.
Jobs1262006/1262007 run native complete optimizations, on the same anode26 with
eight bound workers, BLAS one and unchanged1e-3/1e-7 outer criteria. Saved
1261850/1261856 supply the baseline, without rerunning either executable.

| 240 space | Structures | Steps/HVP, both | Saved SCF / s | New SCF / s | Saved RSS / MiB | New RSS / MiB |
|---|---:|---:|---:|---:|---:|---:|
| COV | 14 | 3/20 | 6.6855 | 6.6624 | 1520.36 | 1523.85 |
| ION(0:1) | 294 | 3/27 | 8.6599 | 8.4470 | 1510.54 | 1516.21 |
| ION(0:2) | 1134 | 4/32 | 10.7693 | 11.8334 | 1547.23 | 1527.26 |
| ION(0:3) | 1694 | 4/32 | 11.7334 | 12.1997 | 1515.35 | 1517.00 |
| FULL | 1764 | 4/30 | 9.7405 | 9.9365 | 1515.57 | 1555.70 |

Maximum240 energy difference is2.0e-12 Eh. The outer trajectory/action counts
are retained, but ION02/03 remain9.9/4.0 percent slower in this native sample.
There is no material whole-process memory reduction. Complete AO integrals,
other accepted arrays and allocator residency remain outside this component.
Do not report a universal speed or peak-memory improvement.

FeCl2 reaches-2181.624464867221 Eh in15 accepted steps/244 HVP,9.7826 s SCF,
415.55 MiB RSS. Saved1261856 is16/281,-2181.624464871859 Eh,12.4227 s,
402.83 MiB. The final energy gap is4.64e-9 Eh. FeCl2 uses the complete-space
inverse, so its trajectory/time change is not evidence of a speedup from
removing restricted completion; no unverified rounding explanation is made.

Separate new-code scalar-phase capture1262012 compares with deletion-only
1261992, not with a rerun of historical code:

| Space | Zero-start M actions | Seeded M actions | Zero-start solve / s | Seeded solve / s |
|---|---:|---:|---:|---:|
| COV | 422 | 333 | 0.05588 | 0.04565 |
| ION(0:1) | 3711 | 2506 | 0.61482 | 0.41325 |
| ION(0:2) | 6148 | 3931 | 1.31769 | 0.85668 |
| ION(0:3) | 7046 | 4371 | 1.66824 | 1.04463 |

All solves retain their true original residual criterion. Instrumented phase
wall is not substituted for native SCF timing. The exact-known-representative
seed removes21--38 percent of selected M actions in these captures, but still
does not recover the older completion's iteration count. Optional completion
was not a useless preconditioner: its conditioning benefit was real, though
its full-space ownership is unnecessary. Further work should express useful
recoupling/conditioning information inside the selected structures, rather
than resurrecting the completed-sector carrier or tuning tolerances.

Machine-readable results: `selected_transport_240.csv`,
`selected_transport_final_fecl2.csv`, `selected_transport_solves.csv`.
`tests/mcscf/summarize_transport.py` parses the raw scalar events. Build and
input/binary/source hashes, affinity, native GNU-time peaks, ctest details and
Slurm accounting are retained with the evidence archive. Initial job1261968
also records a test-fixture failure: topology-only actions cannot construct a
Davidson diagonal without one supplied. The replacement Gram test uses
ActionOnly and explicit tiny exterior minors instead. Corrected fixture and
final production snapshots are qualified by1261985/1262005; no production
runtime error is hidden by this fixture correction.

This eliminates one full-product consumer, not all full-product consumers.
Broad-support accepted/HVP ownership, complete-space inverses and the carrier
AH preconditioner remain independently auditable work.
