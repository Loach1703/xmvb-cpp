# Selected-metric contraction order without closure

2026-10-17. C12/LOFLEA HAO CAS(12,12), exact Libcint/cc-pVDZ,
structure-space Davidson, ISCF=9/CIAH. Native tests use Hanhai22/anode42,
32 physical cores, OpenMP cores/close, BLAS=1 and a 64-GiB Slurm allocation.
No historical executable was rerun. Unrelated local RI edits were not transferred.

## Mathematical and code change retained

For arbitrary selected columns W, the cross metric is
`M_lr v = W^T T_l^T T_r W v`. The intermediate zero extension of `W v`
does not add independent coefficients and does not assume that the selected
space is invariant under T. Both ends still use the original selected W.
In particular, a complete touched product grid is not required to evaluate
this identity by ordered exterior maps.

The old metric and transport kernels admitted this contraction order only
when every product was touched. Other selections used a sum of source
factors, including length-r dots at all t touched products, for each metric
CG action. That O(t r) work is unnecessary when the bounded product-order
contraction is cheaper. The new shared `use_product_overlap()` compares
structural shear/scale/dot work estimates and the simultaneous image payload.
It is used by metric actions, overlap adjoints and the fused transport VJP.
Narrow support and carriers that exceed the existing workspace allowance
continue to use exact factor panels. These are two exact contraction orders,
not different physical models or a numerical fallback.

For this input A=B=r=924, t=242088 and each spin map has 16632 retained
shear edges. The nondirectional structural estimates are 288576288 for
factor order and 64886976 for product order, a 4.447 ratio. This is an
arithmetic estimate for this metric contraction, not a measured kernel
speedup or a prediction for the whole SCF calculation. The common sparse
expansion/projection and hardware-dependent costs are not in this estimate.
For broad balanced support, the avoidable O(t r) dots are replaced by
O(n_active^2 A B) exterior work.

No persistent spin-pair matrix is added. One 924-by-924 double image uses
6.514 MiB; the transport VJP's five images use at most 32.570 MiB under
the unchanged 64-MiB allowance. Images are shared by the thread team, not
replicated per thread. Larger carriers remain panel-streamed. The selected
metric inverse remains an actual selected solve, not a selected block of a
FULL inverse. The existing positive spin-Gram preconditioner, ordinary CI
solver, HVP, orbital coordinates and optimizer strategy are unchanged in
the retained version.

Derivation: `article/source_column_hamiltonian.tex`, Section 10.8,
"Selected-metric contraction order without closure". The PDF was rebuilt
with Tectonic and the added equations and text inspected.

## Native convergence comparison

The three ION(0:2) input files are byte-identical, SHA256
`48b44d64b773d4abaf388fcf728d478a09ab76f5426f889822b1dbc8e031432f`.
Original guesses and ITMAX=2000 are unchanged. Energy tolerance is still
1e-7 Eh, projected gradient tolerance 1e-3, and metric-CG's full true
relative-residual target is still sqrt(machine epsilon).

| Version / job | Space | Structures | Steps / HVP | Final energy / Eh | Application total / s | Peak RSS / GiB |
|---|---|---:|---:|---:|---:|---:|
| Saved FULL / 1265114_2 | FULL | 226512 | 4 / 40 | -399.417495587055 | 102.727418 | 3.572872 |
| Historical ION / 1265093_1 | ION(0:2) | 47256 | 6 / 86 | -399.397232509251 | 1631.498160 | 3.583797 |
| Previous terminal fix / 1265186_1 | ION(0:2) | 47256 | 6 / 86 | -399.397232509252 | 153.175347 | 2.886070 |
| **Retained contraction order / 1265216_1** | **ION(0:2)** | **47256** | **6 / 86** | **-399.397232509246** | **129.879149** | **2.932171** |

The retained change reduces observed application time by **15.21%** versus
the immediately preceding converged version. Both use CPU affinity 96-127
on anode42 and identical inputs. GNU-time process time decreases from
156.06 to 132.79 s. Peak RSS rises by about **47.21 MiB**, consistent with
using bounded shared transient images rather than the old source panels;
it is not a claim of reduced memory. RSS remains below 3 GiB. Measurements
are single native runs on a shared node, not exclusive-node repetitions.

The final energy difference is about +6e-12 Eh. Per-step HVP counts remain
12,12,14,18,27,3. Final projected gradient infinity norm is 4.07731083e-7.
Final CI solve uses 22 H/S block actions in both versions. Thus this change
does **not** claim fewer CI iterations or fewer orbital HVPs. It accelerates
the selected metric/transport contraction schedule. FULL still finishes
faster, with fewer outer steps/HVPs and six final CI actions. Its different
variational energy is not an energy-parity target for ION(0:2).

The last trial has a computed energy increase of 5.22959454e-12 Eh and
passes both original terminal conditions at its own point. It is accepted
without a rejection or added HVP. This native run also exercises the
terminal-trial ordering repaired in the preceding commit.

## Qualification and attribution

Retained qualification **1265213** passes **11/11 tests in 3.47 s**.
`test_representation_overlap` now compares both contraction orders also
for broad but incomplete touched support, using an independent occupied-minor
oracle, directional actions, two-sided adjoints and finite-transport gradients.
Existing tests cover unequal alpha/beta counts, narrow panel budgets and
nonclosed selected spaces. No test tolerance was relaxed.

Frozen bundle **1265214** contains the qualified executable, SHA256
`fd8a01907f317be2336c158c2d345812ba21ddf9376185f94bf19d59ed96b988`.
Native **1265216_1** completes with exit 0. The source manifest and build
configuration accompany the bundle. This is not a later rebuilt executable.

The separate new-code gprofng capture **1265217** is deliberately ITMAX=1
and its application exits 2 at that limit. It is attribution only, not
production timing or a converged run. Collection reports no errors/warnings.
Its first macro uses 12 HVPs; three transport solves need 88,176,89 CG
iterations, with full true relative residuals 1.3903e-8,1.4746e-8,1.4471e-8.
Their measured metric-action wall totals are 0.6531,1.1853,0.6415 s;
the fused transport VJP costs 0.0368 s. No same-input old ION(0:2) phase
capture is available, so no old/new per-kernel speed ratio is asserted.
CI solving, H/S projections and orbital actions remain material costs.
Nested phase scopes are not additive; CPU samples are not wall fractions.

## Experiment not retained

Before isolating the exact contraction change, a candidate also added an
occupation-trace H-E M Davidson denominator and balanced recycling between
the two transport solve batches. Qualification **1265206** passed 11 tests,
but that combination did not qualify as an end-to-end improvement:

- ION(0:1), **1265208_0**: converged, 4 steps / 54 HVP, 42.031273 s,
  -399.323427907220 Eh, 2.672710 GiB.
- ION(0:2), **1265208_1**: cancelled after 13 accepted steps / 133 HVP,
  about 189 s of Slurm elapsed time, without a converged result. The path
  approached the previous endpoint and then continued into a lower-energy
  branch. It is not a completed speedup or energy-parity result.

This combination changed two solver ingredients, so its behavior does not
identify either ingredient alone as the cause. Both additions were removed
from production source, including their unused methods/tests. No optional
switch or compatibility path retains them. Their frozen bundle and diagnostic
logs are kept only as experiment evidence. The retained version keeps both
solvers unchanged and isolates the exact contraction-order change.

## Records and reproduction

`comparison.csv` is generated with `tests/mcscf/c12_spaces.py.summarize`,
with additional input hashes, affinity and final CI counts read from logs.
The parser checks completed convergence and byte-identical ION inputs.
`accounting.txt` records all qualification, freeze, native and capture jobs,
including the cancelled experiment. `qualification.txt` is the retained
11-test result. `evidence.tar.gz` contains current/experimental input bundles
(without executables), native logs, scalar traces, timing/RSS and profiler
reports. Raw `.er` captures and both frozen executables remain on Hanhai22.
Controls reside in the preceding committed evidence archives, not reruns.

The retained live run is
`/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015/c12-spaces-1265216/ion02`.
The new source paths are `structures/assembly/{overlap_tiles,transport_tiles}.cpp`
and the small shared private planner declaration in `action.hpp`, all below
`src/vbscf/`. No new physical-model wrapper or duplicated solver was introduced.
