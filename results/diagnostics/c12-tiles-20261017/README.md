# C12 contraction comparison — Hanhai22, 32 cores

Job **1264541**, `anode40`, completed 2026-10-17, exit 0. Before/after runs use
one allocation, the same input, 32 physical cores, OpenMP `close/cores`, and
single-threaded BLAS. Raw results are in `run-1264541/`; remote binaries remain
under the execution mirror's `c12-tiles-1264541/new/` and the frozen baseline
directory. See `provenance.txt` for binary revisions and hashes. All eight
edited production/test source files match the compute-node source manifest.

## Scope

"C12" here means the existing singlet CAS(12,12), 132-covalent-structure
operator fixture. No separate user-provided `C12.xmi` was found. The real
molecular end-to-end comparison uses LOFLEA's original HAO, exact cc-pVDZ input
and guess, with only `ISCF=7` changed to `9` (CIAH). Numerical thresholds,
workspace budgets and `ITMAX=2000` are unchanged. This is **not** a C14 run.

The before application is the cd958ad7 production executable. The before
kernel executable was last built after 48c145ba; its contractions were unchanged
by the later selected-metric ownership fix. Do not label kernel construction
timing as a strict cd958ad7 measurement. Neither historical source nor the
completed C14 baseline was rebuilt or rerun.

## End-to-end result

| LOFLEA, CAS(12,12), COV 132 | Before | After |
|---|---:|---:|
| Accepted SCF steps | 4 | 4 |
| Orbital HVPs | 36 | 36 |
| Keyframes | 15 | 15 |
| Initial total energy / Eh | -398.533053417516 | -398.533053417516 |
| Final total energy / Eh | -399.106646267332 | -399.106646267331 |
| SCF wall / s | 42.027444 | 41.505578 |
| Application end-to-end / s | 44.996525 | 44.536831 |
| Process peak RSS / MiB, GNU time | 2660.2695 | 2723.3750 |

Both terminated with `ciah_projected_gradient_tolerance`. The measured total
time reduction is **1.02%**, not an established significant speedup; process
peak RSS increased by **63.11 MiB**. Slurm's sampled MaxRSS differs from GNU
time and is retained separately. Batch elapsed time includes compilation,
unit checks, all kernel probes and both molecular runs, not one SCF calculation.

## Fixed-point operator probes

Native in-process phase timings, one cold and one warm action where supported:

| Probe / s | Before | After |
|---|---:|---:|
| H/S action, cold | 0.149600 | 0.206468 |
| H/S action, warm | 0.105439 | 0.127215 |
| Structure response, warm | 0.194061 | 0.189514 |
| Accepted orbital adjoint | 0.305266 | 0.266650 |
| Orbital-direction jet | 0.266877 | 0.285919 |

These short probes have mixed results, including slower H/S action and
orbital direction. They do **not** establish a general operator speedup.
The Gram probe's separately printed `gram_action_seconds` is a **single-RHS**
action and does not isolate the changed multi-RHS selected-metric assembly.

## Numerical checks

- Release compilation passed; **9/9 relevant unit tests passed** in 2.81 s.
- Independent selected H/S and structure-response reference checks passed.
- Structure-response image error: `5.2861e-16`; integral-adjoint error:
  `1.4979e-15` (norms, as printed by the test).
- Product-jet HVP finite-difference error: `5.9108e-11` before and after;
  scalar-gradient finite-difference error: `5.2119e-12` before and after.
- Prepared selected metric versus independently streamed RHSs:
  `5.9763e-16` before, `5.6035e-16` after.
- H/S Rayleigh energy and response checksum agree at the printed precision.

## Interpretation

Source-major multi-RHS contraction removes a genuine redundant transformation
loop without increasing the byte budget or allocating a product-by-RHS array.
The tiny-workspace multi-batch test verifies this reassociation independently.
However, C12 does not reproduce C14's old one-source-column budget collapse;
neither that operation-count argument nor passed numerical tests imply a
measured end-to-end acceleration here. The retained C14 measurement remains
1190.925493 s, 4 steps, 46 HVPs, 6.1037 GiB. No extrapolated new C14 time is
claimed. Subsequent performance work should use these C12 native results as
the immediate comparator and isolate expensive operator phases before further
implementation changes.
