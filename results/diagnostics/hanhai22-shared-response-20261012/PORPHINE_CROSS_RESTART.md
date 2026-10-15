# Fe porphine: cross-program wavefunction restart and CIAH repair

## Scope and settings

All chemistry and builds ran under Slurm on Hanhai22. The comparison uses the
actual CPP triplet endpoint exported in diagnostic `1247963`: 89 inactive
orbitals, CAS(8,11), 469 Cartesian cc-pVDZ AO functions, full spin-triplet VB
space, no ECP, no pseudopotential, and no spatial symmetry. Gradients and
energy thresholds remain 0.001 and 1e-7 Eh. These are bounded two-macrostep
restarts, not converged endpoint comparisons.

The transfer reconstructs individually normalized raw orbitals, projects the
active orbitals out of the inactive space, and uses the same upper Cholesky
active gauge as the direct-CI carrier. The actual carrier is reordered by
explicit occupation masks into PySCF string order. It is not replaced by an
SCF guess or by unrelated CASCI coefficients. A complete virtual complement
is added without changing the core/active wavefunction.

## Independent gradient and representation checks

Job array `1248045` compares PySCF default one-step AH and coupled Newton.
The previous attempt `1248007` was cancelled after a compute-node probe found
only one effective CPU. The corrected script retains Slurm's 32-core binding
but disables OpenMP runtime rebinding, matching the previously verified
PySCF protocol. Both corrected processes report 32 available CPUs and 32
PySCF threads.

- Normalized AO-overlap difference: 1.0547118733938987e-15.
- Transferred full MO Gram error: 7.661228719073631e-13.
- Actual carrier spin-square: 2.000000000000023.
- CPP directional energy derivative: -0.0058411511061091477.
- PySCF directional energy derivative: -0.005841154252947948.
- Absolute derivative difference: 3.1468388e-9.

The raw `CROSS_DIRECTION` log from this job prints the internal `mc1step`
gradient contraction (-0.002920577126473974), which is half the energy
derivative. PySCF `newton_casscf.gen_g_hop` explicitly multiplies this orbital
gradient by two. The archived log must not be silently relabelled or interpreted
as a factor-of-two derivative bug. The current diagnostic script reports both
conventions explicitly.

## PySCF can optimize the lower-energy wavefunction

| Point | Total energy / Eh |
|---|---:|
| Existing PySCF triplet reference endpoint | -2244.782671264980 |
| CPP exported triplet endpoint | -2244.790784197408 |
| PySCF start after wavefunction transfer | -2244.790785879380 |
| PySCF default AH, macro 1 | -2244.790828570880 |
| PySCF default AH, macro 2 | -2244.790831123155 |
| PySCF coupled Newton, macro 1 | -2244.790833291050 |
| PySCF coupled Newton, macro 2 | -2244.790835030458 |

Both PySCF variants retain the lower-energy branch and further lower its
energy while preserving triplet spin. The old reference endpoint is therefore
not a unique attainable energy for these settings. An approximately 8 mEh
gap to that endpoint is not evidence that CPP's energy functional violates the
variational space or that its carrier has spin contamination.

There remains a much smaller 1.681972e-6 Eh cross-provider difference at the
transferred starting wavefunction. Nuclear repulsion agrees to approximately
5e-12 Eh (PySCF 2517.980293225999; CPP 2517.980293226004). The electronic
offset is not yet localized and is not claimed to satisfy 1e-7 Eh agreement.
No integral or optimizer thresholds were adjusted to conceal it.

PySCF's two-step runs are intentionally not converged. Elapsed Python workflow
times are 259.39 s (default AH) and 409.20 s (coupled Newton); these include
the independent gradient construction and repeated PySCF setup. They are not
matched optimizer-only performance timings.

## Repair: retain descent information already computed

The original CIAH restart accumulates an uphill macrostep, although earlier
microstep prefixes have a descending frozen model. Scaling that final uphill
direction cannot repair its nonpositive predicted reduction.

The repair retains the best descending prefix among microsteps left in place
by the existing checks, together with its already computed Hessian image and
last microstep. It replaces the final accumulated step **only when its model
predicts no decrease**. Valid finite negative-curvature steps with positive
predicted decrease are preserved, including those with positive initial slope.
Actual relaxed energy acceptance and its existing thresholds are unchanged.
No HVP, response solve, Hessian approximation, or new stopping condition is
added. Storage is three coupled-coordinate vectors, not replicated integrals
or large CI matrices.

The last AH Ritz residual is not a certificate of a restored prefix; the
returned record is marked approximate when restoration occurs.

## Bounded NEO validation and regression

Job `1248067` is the two-step porphine validation. Both steps restore a
descending prefix and pass the unchanged nonlinear acceptance on their first
trial:

| Macro | HVP actions | Predicted decrease / Eh | Actual decrease / Eh | Energy / Eh | Physical gradient norm |
|---:|---:|---:|---:|---:|---:|
| 0 | — | — | — | -2244.790784197408 | 0.01136091 |
| 1 | 18 | 1.27755032e-5 | 1.42131139e-5 | -2244.790798410522 | 0.00695461 |
| 2 | 19 | 8.45495754e-6 | 8.88453633e-6 | -2244.790807295059 | 0.00595516 |

The first macro still costs 18 HVPs, as in the rejected restart: retaining a
prefix added no HVP. The native program exits 2 because of the intentional
two-step iteration cap, **not** macro energy rejection or convergence. Wall
time is approximately 316 s including setup and final output. This particular
binary used the initial guard that also restored uphill-slope finite steps;
both observed restorations have nonpositive predicted reductions, so the
final narrower guard makes the same choice for these captured steps. The
final source also corrects the diagnostic approximate-step flag explicitly.

The final source was built and checked in regression job `1248081`. F2 and
FeCl2 ran concurrently in separate 32-core Slurm steps; the CIAH unit test
also passed.

| Case | Outer steps | HVP actions | Energy / Eh | Status |
|---|---:|---:|---:|---|
| F2 | 3 | 24 | -198.751155830525 | Converged |
| FeCl2 | 21 | 456 | -2181.624464672223 | Converged |

FeCl2 differs from the existing analytic-diagonal baseline by 4.51e-8 Eh.
No old baseline was rerun. No claim of full porphine convergence or of
singlet cross-program agreement is made: the cancelled current singlet
endpoint was not checkpointed.
