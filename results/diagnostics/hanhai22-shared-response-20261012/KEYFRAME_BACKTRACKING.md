# True-energy CIAH keyframe backtracking

## Evidence and change

Old porphine singlet job 1248745_0 accepted a nonlinear keyframe in attempted
macro 65 with fixed-coefficient energy -2244.739935866974 Eh, below macro
source -2244.739920658597 Eh. Subsequent frozen-model microsteps produced an
uphill trial. The whole-displacement line search ultimately failed, despite
the already evaluated descending point. The problem is loss of that valid
nonlinear point, not justification for relaxing the convergence test.

The macro now retains the lowest true-energy accepted keyframe below its
source. This stores the displacement in the immutable origin chart, its
already accumulated origin HVP image, the last microstep and scalar Ritz
metadata. If the final real-energy trial is rejected, backtrack to this
evaluated keyframe before scaling a ray from the origin. Reoptimized CI energy
must still pass the existing strict descent test; the stored fixed-C energy
is not blindly accepted as a relaxed state energy. A checkpoint is consumed
once, preventing unlimited retries. No extra HVP, gradient evaluation,
convergence condition or tuned threshold is introduced to save it.

This differs from the previously removed quadratic-prefix selection: no old
origin energy model is used to decide whether a prefix is good. The selection
uses a true nonlinear evaluation. At a restored keyframe the model prediction
is anchored to that evaluation, and frozen images/restart directions are
restored consistently. The recorded AH residual is still a last-microstep
quantity, not a certificate for the complete nonlinear displacement.

There is at most one macro-local checkpoint: three coupled-coordinate vectors,
roughly 2 MiB for the triplet. There is no integral/matrix copy per thread,
cross-macro Hessian history or full Hessian construction.

This does not solve every rejection. In the old triplet rejected macro 51,
the evaluated keyframe itself was uphill; it is deliberately not saved as
a descending checkpoint. A different model defect may remain in that case.

## Verification completed

Five independent unit executables passed in compute job 1249059. New tests
check rejection of uphill points, preservation of the lower evaluated point,
matching displacement/image/restart metadata and one-time consumption.

All native cases use 32 bound Hanhai22 cores, exact integrals, original inputs,
the coherent coupled preconditioner and unchanged 1e-3 / 1e-7 tolerances.

| Case | Job | Macros | HVP | Total energy / Eh | GNU wall / s | GNU RSS / MiB |
|---|---|---:|---:|---:|---:|---:|
| FeCl2 | 1249059 | 13 | 257 | -2181.624464862586 | 9.43 | 488.0 |
| F2 | 1249065_0 | 3 | 23 | -198.751155830526 | 0.28 | 17.6 |
| 241 | 1249065_1 | 3 | 21 | -230.720590392836 | 2.65 | 815.2 |
| MnF2 | 1249065_2 | 10 | 137 | -1348.893353223950 | 6.04 | 532.4 |

All converged, exit 0. Macro/HVP counts match the archived retained versions.
None triggers keyframe backtracking, so these prove regression safety, not
the effectiveness of the new path. The FeCl2 energy differs by +2.1147e-8 Eh
from the existing PySCF coupled-Newton reference -2181.624464883733 Eh.
Lower native times/RSS on another node are not claimed as algorithm speedups.

## Pending targeted porphine verification

1249066_0/1 are bounded four-macro singlet/triplet diagnostics on Hanhai22.
They start from the old terminal neighborhoods, reconstructed from the raw
Molden export by reversing our writer's AO permutation only. No AO rescaling,
new basis or ECP is introduced. The frozen executable is from 1249059.

These are **restarts**, not original-guess benchmarks. They do not preserve
the old workspace seed or its rejected microstep sequence, so crossing the
old stopping neighborhood alone cannot prove the backtracking path helped.
Check the `ciah_keyframe_backtrack` event, actual accepted energy, gradient
and subsequent trajectory. Do not present four restarted steps as the total
convergence cost or claim success until the jobs and convergence checks finish.

The original-guess PySCF singlet 1248772_1 remains live (macro 66 at the latest
observation). Both old NEO jobs are terminal, not converged; no original-guess
replacement NEO long run was submitted. Energy equivalence and efficient
porphine convergence remain outstanding.
