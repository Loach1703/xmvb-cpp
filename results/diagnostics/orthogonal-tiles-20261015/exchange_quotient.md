# Spin-exchange quotient: implementation and measurements

## What changed

For identical alpha/beta string spaces and a common exact transpose parity
of the selected VB columns, projected H/S now contracts only one
representative of each exchange orbit in the one-body-generator intermediate.
Off-diagonal representatives carry multiplicity two. Odd-parity diagonal
representatives vanish and are omitted.

Parity is checked once against the actual sparse expansion, with exact
coefficient comparisons. Equal dimensions alone do not activate this reduction.
There is no energy threshold, rank cutoff, selected-space closure assumption,
or memory-based representation switch. Mixed parity and unequal spins still
use all mathematically necessary rows.

The sigma gather returns a representative for the final selected-space
projection, not a full physical CI image. The ordinary sigma interface keeps
its full-image semantics. Independently track overlap panels, since the
quotient Hamiltonian need not request every input column needed for M.

Batch planning reserves only the upper-triangular domain. The workspace
budget stays unchanged; no per-thread carrier or new full product array
is allocated. This is a roughly twofold reduction of one intermediate,
**not a new asymptotic complexity and not a complete HVP optimization**.

Theory: [spin_exchange_contraction.md](../../../article/spin_exchange_contraction.md).
Separate integral-adjoint/response kernels have not yet been quotient-reduced.

## Verification and provenance

- Parent commit: ea211130.
- Hanhai22 only; gcc14.1.0/cmake3.30.5, explicit core binding, BLAS1.
- 240 final job1260261, anode15, 32cores; all six CTests and five
  converged cases passed. Initial candidate job1260255 is retained in CSV.
- Added mixed-parity equal-dimension test compiled and passed in
  job1260276, anode26. The suite includes +/- parity, multi-RHS, nonclosed
  selected spaces, unequal spins, and narrow streaming budgets.
- CAS12 job1260276 uses8cores on anode26, matched to saved job1259764
  (8cores verified via Slurm accounting). Historical binaries were not rerun.
- Complete H/S images agree with the independent reference to
  2.14e-14 absolute H error and 2.02e-15 M error in the CAS12 cases.
- Raw native output, time-v RSS, binary/source hashes and unit-test output
  are retained with each job. Native and reference-process timings are
  recorded separately by summarize_selected_tiles.py.

## CAS(12,12) complete H/S action, 8 cores, same node

This is a fixed-orbital synthetic kernel, not a molecular SCF optimization.
Native warm action time is shown. RSS is whole-process GNU-time/getrusage
peak, not exclusively the sigma workspace.

| Selected structures | Previous ms | New ms | Change | Previous RSS MiB | New RSS MiB |
|---:|---:|---:|---:|---:|---:|
| 1 | 87.31 | 59.80 | -31.5% | 63.59 | 68.28 |
| 3 | 87.78 | 72.84 | -17.0% | 65.66 | 68.27 |
| 20 | 102.60 | 90.93 | -11.4% | 66.69 | 70.39 |
| 132 | 100.11 | 113.49 | +13.4% | 67.69 | 71.66 |

The last result is not hidden. The separate new correctness-reference process
measured88.68ms before constructing its independent reference. A subsequent
same-node, same8core stage-timing run1260278 measured88.60ms. Thus the
132-structure speedup/regression is not robustly established by these short,
nonexclusive-node timings. Do not substitute the fastest value into the
native comparison table.

Stage timing for132structures confirms304590 generator representatives
instead of608751 original closure points:429 fixed points and two-member
off-diagonal orbits. Batch count falls from32 to17, source columns from11722
to8081, and source elements from10642550 to6908066. These are work counts,
not a cross-core timing comparison (saved stage instrumentation used32cores).
Peak reported batch payload increases from34.90MB to37.45MB, still within the
unchanged64MiB target. More beta groups share each batch; RSS has not improved.

The new8core warm phase totals are prepare8.55ms, read2.21ms, compute31.07ms,
write28.12ms. Input preparation and final projection also remain outside these
sigma subphases. Removing half the generator rows does not halve these other
costs. A gprofng availability probe on the compute node failed; no sampling
percentages or function-level hotspot claims are made.

## 240 convergence, 32 cores, same anode15

Compared with the immediately preceding source-factor version1260223.
Exact/libcint, cc-pVDZ, NEO/Davidson, analytic diagonal, gradient1e-3 and
energy-change1e-7Eh are unchanged. Both methods have identical steps/HVP.

| Space | Steps/HVP | SCF seconds previous→new | Process seconds previous→new | RSS MiB previous→new | Final energy Eh |
|---|---:|---:|---:|---:|---:|
| COV | 3/20 | 2.133→2.061 | 3.64→3.59 | 1546.4→1539.5 | -343.239211501261 |
| ION(0:1) | 3/27 | 3.078→3.082 | 4.65→4.56 | 1607.3→1531.3 | -343.407736636061 |
| ION(0:2) | 4/32 | 4.041→4.032 | 5.62→5.58 | 1544.6→1597.1 | -343.443415326851 |
| ION(0:3) | 4/32 | 4.157→4.188 | 5.65→5.77 | 1566.9→1527.3 | -343.446256564846 |
| FULL | 4/30 | 3.391→3.388 | 4.94→4.89 | 1613.6→1534.6 | -343.446302545201 |

Final energies agree within2e-12Eh. End-to-end timing is effectively flat;
this is not evidence of a large molecular speedup. Full generator target
count changes from4900 to2485, but other optimization costs remain.

## Pending

At report creation, the32core CAS12 kernel job1260267 and dependent stage/
profiler job1260268 are queued for anode15 (the saved32core stage baseline
also used anode15). They use the same final numerical source. No successful
32core kernel timing or sampling result is claimed yet.

Next work should target the still-unreduced integral adjoint/response and
projection chains, retaining this exact parity identity. The current result
is a verified smaller contraction and a few-structure kernel speedup, not
universal speedup or solved large-active-space scaling.
