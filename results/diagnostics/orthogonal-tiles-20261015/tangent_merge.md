# Exact merging of orbital and coefficient tangents

Parent: `920623fc`. Builds and numerical work run only in Hanhai22 Slurm,
anode26, 8 bound OpenMP cores and single-thread BLAS. Compilation flags and
thresholds match the previous checkpoint.

## Algebra and implementation

Each coordinate-anchored source has one constant raw spin factor. Its
coefficient derivative can therefore be added to the orbital derivative of
the opposite factor. `SpinSources::tangent_source` represents the combined
direction with two factor families rather than concatenating a third family.
It uses one single-spin scratch panel, overwriting the generator image after
accumulation. Exchange-parity reconstruction changes anchors and signs
consistently. No SVD, amplitude cutoff, closure assumption or optimizer policy
change is introduced. Product-order arithmetic estimates use the new exact
factor count as well.

Equations52--55 of `article/source_column_hamiltonian.tex` document the proof.
The rebuilt PDF was inspected; there are no overfull or undefined-reference
warnings.

## Verification

- Job1261360:9/9 focused tests pass. New independent full-carrier generator
  reconstruction covers mixed anchors, unequal spins, both parity signs,
  coefficient-only/orbital-only/combined directions and panels crossing every
  family boundary. Actual12/14-orbital source HVP consumers agree with the
  separate complete-carrier reference; their full values are in
  `tangent_state.csv`.
- Job1261367:9/9 focused tests and9 complete240 SCF regressions pass, including
  arbitrary1/3-structure selections with dense and Davidson CI. Steps and HVP
  counts remain unchanged. Maximum energy difference from the saved baseline
  is3.98e-12 Eh. See `tangent_240.csv` and `tangent_selected_240.csv`.
- Job1261368: actual production-gradient, orbital/coupling-HVP and central
  relaxed-gradient-difference probe on the constructed exact16-AO CAS16/16
  Hamiltonian. This is not a molecular convergence benchmark.

|Production probe|Saved1261340|Merged1261368|
|---|---:|---:|
|Gradient / s|25.731710|25.898488|
|HVP / s|153.490708|154.537731|
|Full probe, including difference points / s|230.011520|231.269835|
|Peak RSS / MiB|155.6953|155.9414|
|HVP difference error|1.7278e-9|1.7278e-9|
|Accepted raw/orthogonal/residual product doubles|0/0/0|0/0/0|

The merged formula is exact but **has not produced an end-to-end speedup**.
It removes duplicated factor families, not all repeated generation of their
values. No speedup is inferred from the factor count alone.

The production HVP now reports internal stage times:

|Stage|Seconds|
|---|---:|
|Active integral direction|0.000133|
|Selected structure direction|55.008766|
|Active backward|99.527373|
|Structure-response orbital pullback|0.000335|
|Raw product response reconstruction|0|

These are direct stage timers, not sampled CPU percentages. General nested
timers must not be added indiscriminately.

## Attribution tools

Job1261369 stopped before a workload capture because gprofng is absent.
Job1261385 used the installed VTune software sampler, but collection was denied
by the node ptrace policy; its perf probe also failed. No security settings
were changed. Neither attempt supplies a CPU hotspot profile. The failures
are preserved, and an opt-in current-code phase/work-counter capture1261387
is running separately. Native timings above do not include those capture
attempts.

Next examine repeated directional-source generation and the placement of the
physical parity projector in the symmetric integral adjoint. The known
within-state normalization response may also admit exact reuse of accepted
adjoints. These are algebraic opportunities to validate, not measured gains.
