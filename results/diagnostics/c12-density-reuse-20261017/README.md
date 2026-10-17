# Validated density identity, rejected full production rewrite

The independent identity is recorded in
[normalization density reuse](../../../article/normalization_density_reuse.md).
The one-electron integral cotangent is the symmetrized one-body generator
moment, also for distinct nonstationary carriers and unequal spin sectors.
The new independent test is retained in `test_orthogonal_direct_ci`.

The complete candidate reused that cotangent in normalization terms, stored
Y/HY in the same compact two-component cache and omitted unused reads. It
removed 64 generator calls (704 -> 640), reducing that scoped timer from
0.380830s to 0.272936s. HVP finite-difference error remained 5.84911563e-11;
all 9 regression tests passed in 3.45s. However the complete native runs did
not improve. Therefore none of the production cache/consumer changes remain.

Hanhai22, anode40, identical physical cores60-91, original byte-identical
LOFLEA CAS12/COV132/HAO/exact cc-pVDZ, ISCF9 and Davidson. The second native
run replayed only the frozen NEW candidate, without compilation or an old run.

| Metric | Saved current 40595b4e | Candidate 1264640 | New-candidate repeat 1264642 |
|---|---:|---:|---:|
| Final energy / Eh | -399.106646267331 | -399.106646267331 | -399.106646267331 |
| Steps / HVP | 4 / 36 | 4 / 36 | 4 / 36 |
| SCF wall / s | 28.207024 | 29.361272 | 28.839069 |
| Application wall / s | 30.863387 | 32.108210 | 31.809545 |
| Peak RSS / KiB | 2892468 | 2721076 | 2706596 |

Separate phase runs are 31.477191s saved versus 32.582592s candidate. Physical
sigma/adjoint was not faster (10.338633s -> 10.762344s), despite fewer generator
calls. Unchanged scopes also vary; no unmeasured cause is asserted. A faster
small contraction is not treated as sufficient evidence for faster convergence.

After removing the candidate, job **1264646** rebuilt the canonical production
code plus the new independent density-identity test on a Hanhai22 compute node
and reran **9/9 tests**, all passing in **2.70s**. This was a unit/build check,
not another SCF timing or an old-version benchmark. Source hashes are retained.
The canonical cache continues to hold Y/(H-E)Y; no runtime switch, legacy
consumer, added convergence check or weaker tolerance remains.

The main remaining direction is a dataflow change rather than another layout
tweak: [two-slab directional sigma](../../../article/two_slab_direction.md).
It proposes reusing the dead halo as shared GEMM micro-panel scratch, then
sharing an accepted scatter across orbital and CI requests in temporal order.
This is NOT implemented or qualified yet; it must keep separate coarse images
and the existing global memory budget. The draft states the required tests.
