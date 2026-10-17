# C12 projection review: minimal reuse retained, other prototypes removed

Parent:55b5b17f. All new jobs ran on Hanhai22 anode40 through Slurm,
32 physical cores, affinity60-91, OpenMP32/BLAS1. C++23 Release flags,
dependencies, input, basis, original guess, thresholds and ITMAX2000 are
unchanged. The launcher is tests/mcscf/c12_regression.sbatch. No old binary
was rebuilt or rerun;1264786 is saved evidence.

Scope: the existing C12 operator fixtures plus converged **LOFLEA CAS(12,12)**,
COV132/HAO/exact Libcint/cc-pVDZ, Ualpha=Ubeta924, nBF150, CIAH/Davidson.
These timings are not a separate carbon-C12 molecule or the earlier C14 input.
Input SHA256:092d545498678b0cf772c5c84b9eab674b694375558cd391a1b739b7bdcc35eb.
Basis SHA256:d53e7e02d55fb99cdf823d50cf13fec375e0cbe1220ff61816aee7ab38e38ec9.

## Native complete convergence

All rows converged in4 accepted steps,36 orbital HVPs and36 physical structure
responses. No physical response or residual/gauge term was screened away.

| Candidate / job | SCF / s | Application / s | Peak RSS / MiB | Final energy / Eh |
| --- | ---: | ---: | ---: | ---: |
| Saved parent /1264786 | 21.394834 | 24.205419 | 2650.95 | -399.106646267332 |
| Combined product copyback /1264789 | 21.659273 | 24.503898 | 2697.60 | -399.106646267331 |
| Borrowed projection/workspace /1264791 | 22.183739 | 24.900395 | 2685.46 | -399.106646267331 |
| Source-row packing /1264792 | 20.570733 | 23.404603 | 2702.48 | -399.106646267330 |
| Supported projection jet /1264796 | 21.264987 | 24.068527 | 2694.72 | -399.106646267330 |
| Supported jet + accepted density /1264804 | 21.664261 | 24.767611 | 2651.13 | -399.106646267331 |
| Minimal accepted density /1264805 | 21.421274 | 24.221871 | 2761.80 | -399.106646267331 |

Process RSS is GNU-time KiB/1024, not the sum of per-thread memory or the Slurm
batch-shell maximum. Native final energy difference is1e-12Eh for the retained
candidate; HVP counts per step remain6,17,10,3. End-to-end time is effectively
unchanged. This round establishes **no significant native speedup or RSS
reduction**. The measured RSS is higher despite adding no numerical workspace;
its cause was not profiled and is not assigned to allocator/thread noise.

## What remains in production

Only accepted-state normalization density reuse remains. The one-electron
integral cotangent already is D(Y,Y); apply its -2wE multiplier to the upper
chart and reuse it for the single-state gauge density. Remove the separate
accepted generator contraction and omit a later carrier reread when neither
images nor diagnostics need it. Oriented nonstationary transition moments
are not replaced by a symmetric density. The Y/(H-E)Y cache, its budget,
all directional contractions and the solver are unchanged.

This is deliberately narrower than the earlier rejected full Y/HY-cache
rewrite in c12-density-reuse-20261017. That old result was not rerun.
The independent mathematical identity is in article/normalization_density_reuse.md.
The other retained change is a no-independent-forcing temporal regression
test against non-tiled image/integral-adjoint oracles.

## Rejected paths and separate attribution

- Combining two necessary Hamiltonian products in shared scratch removed a
  copyback, not a product; complete coupled time stayed7.98s versus7.99s.
- Borrowing blocks and reusing projection allocations did not improve native
  runtime. No runtime switch or alternate allocation path remains.
- Source-row packing is an exact post-transform read-support restriction,
  not a selected-subspace closure assumption. However projection time grew
  from1.828959479s to2.484453354s. Its native23.40s alone is not sufficient:
  the separate instrumented complete run was24.631926s versus24.682620s.
- The directional support jet retained zero-primal/nonzero-direction shear
  edges and passed unequal-spin/nonclosed-seed tests. Nonetheless directional
  transform time was0.749981353s versus0.739086715s and total time did not
  improve. Its additional API, validation refactor and production branch were
  removed. Prototype source/tests are archival evidence only.

Final1264805 separate phase run: application23.292908s, coupled images7.630408856s,
projection1.769096859s. Generator scope is352calls/0.253964164s, versus the
saved416calls/0.312359259s. The64 eliminated calls are a deterministic work
reduction; the approximately0.058s scoped saving cannot justify claiming that
the larger whole-run difference was caused by density reuse. Scope totals
overlap and are not additive time shares. Native measurements above are the
production comparison, not instrumented runs or operator microbenchmarks.

## Validation and retained evidence

Final job1264805 completed;9/9 regression targets passed in3.60s. The extra
CAS12 parity, structure-response reference, product-jet finite-difference and
temporal fixtures also passed. All final production/test hashes match the
recorded compute build. The newly introduced non-retained transform files
were restored before this final build. Existing tests cover nonclosed selected
support, unequal spins, exchange parity and nonstationary response adjoints.

comparison.csv is generated with the existing saved-run row() parser in
c12-temporal-responses-20261017/summarize.py. phase-summary.json uses
tests/mcscf/summarize_c12_cost.py. Both only parse saved data.
evidence.tar.gz contains all six new jobs' inputs, full outputs, convergence
traces, operator checks, time/RSS, build flags, affinity, source/binary hashes
and rejected candidate source/math notes; frozen binaries are excluded and
remain on Hanhai22. The three pre-existing unrelated dirty files were neither
transferred nor staged. No new Ualpha*Ubeta array or per-worker carrier exists
in the retained change. The remaining major cost is still the physical channel
Hamiltonian products and selected output contractions, not the removed norm
moment. Those require a larger demonstrated dataflow improvement.
