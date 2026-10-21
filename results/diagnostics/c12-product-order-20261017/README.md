# Ordered transforms in the executed product-tile chain

2026-10-17; parent `b117bb4d`. Exact Libcint, cc-pVDZ, HAO CAS(12,12),
ISCF=9/CIAH, structure-space Davidson. Neither optimizer policy nor any
energy, gradient, linear-solver or numerical screening threshold changed.
Unrelated local RI changes were not transferred or included.

## Mathematical change and ownership

The preceding source-amplitude changes did not accelerate these broad-support
examples because they execute `contract_product_tiles` and `product_jet`.
This round changes those actual consumers.

For arbitrary selected structures, zero-scatter the original support into Z.
The exact forward map is Y = Ta Z Tb^T and its reverse is
W^T vec_Omega(Ta^T F Tb). Closure of the selected structure space is not
assumed. The structure-space independent coefficients are unchanged.

The old beta-row panel order repeatedly multiplied/sampled every original
product. Across reachable panel columns its beta work is O(G * Ub_sampled).
An ordered exterior action instead costs O(Ua * (Ub + Lb)), where Lb counts
beta shear edges. `prefer_ordered_beta` compares these arithmetic estimates,
counting only reached panel columns. It is not a wall-time oracle. COV's
narrow support keeps row sampling; broad ION/FULL can use ordered factors.

H/S actions now expand once per RHS batch, transform once, let sigma halos
borrow the same columns, accumulate the image, and pull back once. These
input/output owners replace the numerical LRU payloads inside the same
workspace quarters. Sigma no longer needs another overlap carrier: its
finished input becomes the overlap pullback. No complete exterior matrix,
global numerical U-squared cache, or per-thread carrier copy was added.

The accepted/directional product jet can use the same ordered beta contraction
on broad but incomplete raw support, not only on a complete raw product grid.
Both exterior derivative terms and the selected-space nonstationary residual
remain. Sparse accepted-image retention gathers directly into its destination
using exact offsets/rows; there is no temporary packed panel. For carriers
that exceed the existing shared residency budget, the exact streamed order
remains available. This is a residency/contraction-order choice, never an
admission test for orthogonal CI. It does not reduce every large-carrier
streamed case to the resident arithmetic bound.

The derivation and complexity are in Section 12 of
`article/source_column_hamiltonian.tex`; its 33-page PDF was rebuilt and
the new formula pages inspected.

## Correctness qualification

Final build/test job **1265582**, anode06: **14/14 passed**, 6.64 s.
The tests include independent H/S images and finite differences, multi-RHS,
tiny-budget streamed order, unequal spins, both exchange parities, incomplete
broad support with interior column holes, actual compact accepted-image
retention, orbital directions, joint orbital/CI response and generator jets.
The new broad-support fixture explicitly requires the product path, rather
than allowing a passing source-order comparison to hide routing mistakes.

First development qualification 1265565 exposed a bug in a new test fixture,
not an operator crash: it cleared a cache and called an accessor that did not
rebuild it. Temporary stack capture 1265568 located that test line. The fixture
now calls the real accepted-state construction entry. Qualification 1265573
passed before independent review removed the extra gather temporary, refined
the work estimate and added exchange-parity coverage. Superseded current-code
array/profile 1265575/1265576 were cancelled; they are not benchmark results.
No historical executable was rerun.

## End-to-end convergence: saved excellent controls, not old reruns

Final frozen bundle **1265583**, native array **1265585_[1,2,6]**, Hanhai22,
32 physical cores, OpenMP cores/close, BLAS1, 64 GiB requested. All three
inputs match the saved controls byte for byte, including the original guess.
No new SCF iteration/time cap was imposed. Each run completed with exit 0
and `ciah_projected_gradient_tolerance`; there are no failed-HVP actions.

Application time uses the historical CSV definition: initialization plus SCF.
Process time additionally includes final output/teardown, important for FULL's
large structure tables. RSS is GNU time's peak process RSS, not summed Slurm RSS.

| Space | Structures | Steps old/new | HVP old/new | Application s old/new | Process s old/new | RSS GiB old/new |
|---|---:|---:|---:|---:|---:|---:|
| COV | 132 | 4 / 4 | 36 / 36 | 17.208 / 17.838 | 17.40 / 18.06 | 2.572 / 2.599 |
| ION(0:2) | 47,256 | 4 / 4 | 44 / 44 | 104.286 / 89.414 | 107.11 / 92.24 | 2.824 / 2.877 |
| FULL | 226,512 | 4 / 4 | 40 / 40 | 97.037 / 94.786 | 107.31 / 105.26 | 3.568 / 3.522 |

| Space | Current energy / Eh | Current minus saved / Eh |
|---|---:|---:|
| COV | -399.106646267331 | +2.05e-12 |
| ION(0:2) | -399.397232509244 | +9.66e-13 |
| FULL | -399.417495587056 | -1.99e-12 |

Controls are COV **1264892**, ION **1265266_1**, FULL **1265093_2**.
For FULL this deliberately uses the faster saved **97.037-s** result, not
the later 102.727-s control. All nodes are AMD EPYC 9654. Affinities differ:
COV anode40/0-31 -> anode07/96-127; ION anode65/96-127 -> anode06/0-31;
FULL anode42/32-63 -> anode06/32-63. This is one native run per space,
not repeated statistics or a same-allocation causal experiment.

ION saves **14.26%** application time against its excellent control, and
**25.93%** against the immediate parent's saved 120.715 s. FULL saves 2.32%
against its best saved control; this small difference is not asserted as a
robust speedup. COV remains 3.66% slower than its excellent control, despite
being faster than the immediate parent's 19.475 s. No acceleration of COV's
unchanged narrow beta contraction is claimed. There is not an all-space win.

## Actual-path attribution, not another native timing

Current short capture **1265586**, anode07, uses gprofng and phase timing on
one ION macrostep. It intentionally returns 2 at that diagnostic step limit;
the three native convergence runs above all return 0. Capture reports no
errors/warnings. Existing capture **1265545**, anode04, supplies the previous
comparison; it was not rerun. Nested phase walls are not additive and profiler
percentages are not treated as precise useful-compute fractions.

| Instrumented phase | Calls old/new | Wall s old/new |
|---|---:|---:|
| Product generation | 692 / 196 | 4.683 / 1.245 |
| Product projection | 1390 / 317 | 4.104 / 1.621 |
| Left exterior transform | 3223 / 1792 | 1.953 / 0.937 |
| Right exterior transform | 1215 / 1697 | 1.708 / 1.967 |
| CI sigma | 151 / 151 | 8.961 / 6.472 |
| Orbital HVP | 9 / 9 | 3.156 / 2.237 |
| CIAH macrostep | 1 / 1 | 20.502 / 14.130 |

Logs explicitly show `XMVB_PRODUCT_ACTION,...ordered_beta=1`, a single
Hamiltonian projection per batch and zero input/output evictions. The right
ordered transform count increases: that work replaces much more expensive
per-panel original-support reconstruction and sampling. HVP and sigma counts
are unchanged. These are executed product-chain changes, not an inference
from a source kernel that the native jobs never call.

## Reproduction and evidence

Run `python3 results/diagnostics/c12-product-order-20261017/summarize.py`
after unpacking `evidence.tar.gz` here to regenerate `native.csv`,
`comparison.csv` and `phases.csv` from outputs and saved controls. The script
checks convergence, failed actions and all three input hashes before comparing.
The archive retains qualified source/build manifests, exact inputs, complete
native outputs/traces/resource records and profiler reports, excluding binaries
and the large capture directory. The original current capture remains at
`/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015/c12-ionic-profile-1265586/capture.er`.
`accounting.txt` retains Slurm states, including development/cancelled jobs.
