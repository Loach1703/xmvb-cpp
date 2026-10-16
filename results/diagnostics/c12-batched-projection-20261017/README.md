# Batched selected-structure projection and accepted left jets

Parent:20c9e258. Production changes are confined to
src/vbscf/structures/assembly/product_jet.cpp; regression extensions are in
tests/vbscf/unit/test_structure_tiles.cpp. The derivation is recorded in
article/batched_structure_projection.md. Three unrelated user-dirty files
were neither transferred nor staged.

## Strictly shared work

1. Independent orbital/CI Hamiltonian images share one alpha transform call,
   beta row resolution and selected-product traversal. Outputs remain separate.
2. Load both overlap directions directly into adjacent shared panels. Finish
   all oriented generator consumers before their destructive projection.
3. Unpack accepted [Y,(H-E)Y] once into a disposable panel, reconstruct its
   second component in place and project the left jets of [Y,HY] together.
   Remove a residual copy, a separate HY temporary and the later dense Y read.
4. A by-value projection consumes final-use owners and avoids the former
   block-expression temporary followed by another internal owning copy.
   Overlap capture/norm consumers run before final physical projection,
   avoiding another copy of both directions. Accepted cache owners are never
   moved or mutated.
5. Broad resident outputs without a left-direction component use an ordinary
   beta adjoint. Do not allocate/transform a zero carrier for CI-only or
   integral-forcing outputs. This is an operator-component dependency, not
   numerical screening or a memory-based wave-function branch.

The selected subspace need not be closed. No Hessian component, overlap gauge,
nonzero CI-residual term or structure response is dropped. Panel geometry,
accepted Y/(H-E)Y residency, solver policy and accuracy are unchanged. No
per-thread carrier or new full-product array is introduced. Two image blocks
reuse the existing leading workspace allowance, rather than enlarging a tile
or adding a cache. The formal Hamiltonian arithmetic order is unchanged.

## Matched native complete convergence

Hanhai22 anode40,32 physical cores/affinity60-91, OpenMP32/BLAS1, same GCC14.1.0
and CMake3.30.5, C++23 Release -O3 -DNDEBUG -g1 -fno-omit-frame-pointer,
Eigen BLAS/LAPACKE. Existing C12 operator fixtures accompany full original-guess
LOFLEA CAS(12,12), COV132/HAO/exact Libcint/cc-pVDZ, Ualpha=Ubeta924, nBF150,
CIAH/Davidson. This is not a separate carbon-C12 molecule or the C14 input.
Original thresholds and ITMAX2000 are retained. Input and basis hashes match
the saved runs; raw input hashes are archived with each job.

| Version / job | Steps | HVP | Final energy / Eh | SCF / s | Application / s | Peak RSS / MiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Saved row-Gram /1264786 | 4 | 36 | -399.106646267332 | 21.394834 | 24.205419 | 2650.95 |
| Saved immediate parent /1264805 | 4 | 36 | -399.106646267331 | 21.421274 | 24.221871 | 2761.80 |
| Output batch /1264813 | 4 | 36 | -399.106646267332 | 20.279518 | 23.174799 | 2708.91 |
| Plus accepted left-jet batch /1264815 | 4 | 36 | -399.106646267329 | 19.777675 | 22.754901 | 2698.57 |
| Core fusion, omit absent left component /1264816 | 4 | 36 | -399.106646267329 | 19.622553 | 22.417967 | 2635.41 |
| Same1264816 binary, native check /1264817 | 4 | 36 | -399.106646267330 | 20.381447 | 23.075239 | 2757.49 |
| Final capture-consumer ordering /1264818 | 4 | 36 | -399.106646267330 | 21.020280 | 23.819179 | 2704.81 |

The core-fusion measurements are4.7--7.4% shorter than the saved immediate
parent. The last version reorders capture consumers for lower diagnostic
scratch needs; its native application is only1.7% shorter than the parent.
Do not quote the fastest candidate as the final version or as a universal
speedup. No cause is assigned for the between-run timing spread. RSS is GNU
time process maximum, KiB/1024. Its range overlaps the saved measurements;
no stable process-RSS reduction is claimed. There is no extra numerical
workspace in the retained algorithm. Four steps and36 physical structure
responses are unchanged; per-step HVP/response counts are6,17,10,3.

The native check was requested because instrumented full-run time varied.
It executed exactly the then-current NEW1264816 binary, verified by its
executed SHA256, before the final diagnostic-consumer edit. Final1264818 was
separately rebuilt/tested and run to convergence. No old version was run.
All native runs are
full convergence, not one-step extrapolations or operator-only timings.

## Separate phase attribution

| Scope/counter | Saved parent1264805 | Output1264813 | Cached jet1264815 | Core1264816 | Final1264818 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Application / s | 23.292908 | 24.504082 | 23.416355 | 24.415606 | 24.085821 |
| Product projection calls | 1500 | 1212 | 1068 | 1068 | 1068 |
| Product projection / s | 1.769096859 | 1.811163002 | 1.463913893 | 1.504778482 | 1.496966856 |
| Ordinary left-transform calls | 3325 | 3037 | 3037 | 3037 | 3037 |
| Directional-transform calls | 760 | 760 | 616 | 616 | 616 |
| Complete coupled images / s | 7.630408856 | 7.918041416 | 7.351532363 | 7.635223200 | 7.599648872 |

Projection scope is15--17% shorter in the completed fused-jet candidates,
with432 fewer projection passes. Logical independent outputs are still all
produced; this is batching, not fewer HVPs. The final complete coupled timer
is effectively unchanged from the immediate parent, despite the smaller
projection scope. Do not assign the full native gain to that one scope or
claim the instrumented application improved: it did not. Scopes overlap and
instrumented timings are not the native production comparison.

All candidates keep129 generated sigma sweeps,1171 channel batches and
1,120,445,436 reader elements in those sweeps. The eliminated later Y read is
outside that counter. The remaining channel multiplication, scatter and
output work is still present; no asymptotic exponent reduction is claimed.

## Qualification and evidence

Final job1264818 passed9/9 regression targets in3.21s. Added joint tests cover full
touched support with resident projection, budget128 eviction/replay, accepted
cached and uncached images, unequal spins, two roots, and capture on/off.
Existing tests include both exchange parities, zero directions/weights,
nonclosed selected spaces and independent complete-carrier response oracles.
The extra selected12, structure-response reference, product-jet finite-
difference and temporal12 fixtures also passed. Representative errors:

- structure response image3.40e-16; integral adjoint1.27e-15;
- product-jet HVP finite-difference2.5934e-11;
- Lagrangian directional derivative6.9466e-12.

These tests establish operator equivalence at their existing tolerances,
not bitwise equality of every intermediate optimizer coordinate. Native
final energies differ from the saved parent by at most2e-12Eh.

All final compiled production/test hashes match the local retained files.
evidence.tar.gz contains all five new jobs' raw input/output, CMake flags,
tests, traces, time/RSS, affinity, source/binary hashes, the two staged candidate
sources and the final source/test snapshot, excluding frozen binaries.
Final binaries remain under c12-regression-1264818/bin on Hanhai22. sacct.txt records
terminal scheduler states. native-repeat.sbatch records the one current-only
confirmation. comparison.csv uses the existing saved-run row() parser;
phase-summary.json uses tests/mcscf/summarize_c12_cost.py. No old code was rerun.
