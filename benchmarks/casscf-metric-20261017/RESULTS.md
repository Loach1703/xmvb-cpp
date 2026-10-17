# Isolated AH span-metric experiment

Parent: `5c2ad258`. Hanhai25 Slurm, 32 physical cores, exact integrals,
single state. No old binary or PySCF reference was rerun. Build 291978
passed all 12 selected tests with the additional metric assertions.
Runs 291979–291981 all converged. Timings below are single native samples.

## Derivation and controlled change

The VB inactive/active span metric reduces to weights 2 for inactive–active
rotations and 1 for the other two independent orbital sectors in the
orthogonal complete-CAS chart. The existing sqrt(2)-scaled CI tangent has
unit metric. See `article/casscf_ah_metric.md` for the projector derivation.

The experiment changed only the AH metric, matching coarse metric images,
the metric-dependent shift in the existing diagonal preconditioner and the
existing keyframe length callback. Gradient/HVP, coordinate progress norm,
orbital/CI finite updates, seeds, coarse algorithm, microstep decisions,
component clipping and all convergence thresholds remained unchanged.
The .03 component cap does **not** gain an additional sqrt(2) factor: the
VB and CASSCF physical rotation component caps already agree in this limit.

## FeCl2 native `.xmi`

CAS(10e,8o), quintet, 136 Cartesian def2-TZVP AOs, 25 inactive orbitals.
The same input bytes and AO-loader convention were used in both runs;
both ran on anode029 (AMD EPYC 9755). Initial energy is
-2181.617631912270 Eh.

| AH metric | Steps | HVP | Final energy / Eh | Solver / s | End-to-end / s | Peak RSS / MiB |
|---|---:|---:|---:|---:|---:|---:|
| Identity, stored 291946 | 18 | 255 | -2181.624464862390 | 3.502806 | 3.794071 | 407.47 |
| Span, 291981 | 20 | 267 | -2181.624464852350 | 3.839495 | 4.146485 | 402.25 |

The energies differ by 1.00e-8 Eh. Solver wall increases 9.6%; the actual
outer/HVP counts also increase, so this cannot be dismissed as timing noise.
The few-MiB peak difference does not establish an algorithmic memory gain.

## FeCl2 identical numerical AO/MO fixture

Both runs read the unscreened fixture with SHA-256
`801c0fddf648e281e0453fe00a0813a525514b98ff047908fbaad36594232190`.
Both ran on anode034 (AMD EPYC 9755). Initial C++ energy is
-2181.6176319123533 Eh. This path uses the larger unscreened AO graph and
must not be mixed with native-loader wall/RSS.

| AH metric | Steps | HVP | Final energy / Eh | Solver / s | Process / s | Peak RSS / MiB |
|---|---:|---:|---:|---:|---:|---:|
| Identity, stored 291944 | 16 | 209 | -2181.624464858833 | 12.105322 | 12.81 | 1664.00 |
| Span, 291979 | 14 | 204 | -2181.624464854085 | 11.200599 | 11.96 | 1664.00 |

Here the metric helps: two fewer outer steps, five fewer HVPs and 7.5%
lower solver wall, with a 4.75e-9 Eh final-energy difference. It does not
establish a representation-independent improvement.

## F2 small-system regression

| AH metric | Steps | HVP | Final energy / Eh | Solver / s | Process / s | Peak RSS / MiB |
|---|---:|---:|---:|---:|---:|---:|
| Identity, stored 291953 | 3 | 14 | -198.761111541649 | .143675 | .17 | 56.00 |
| Span, 291980 | 4 | 18 | -198.761111541650 | .214160 | .24 | 53.92 |

The final energies agree, but the iteration and action counts regress.
The nodes differ for this small-system test, so timing alone is not used
to qualify that regression.

## Decision and interpretation

The span metric is mathematically well defined and cheap, but this isolated
factor does not explain the faster VB trajectory or consistently improve
CASSCF. It changes the finite AH path; it is not equivalent to transferring
VB's complete horizontal lift, CI transport and nonlinear model. In
$(H-\omega G)d=-g$, the explicit influence of $G$ also vanishes as
$\omega\to0$ in the Newton limit; a different metric cannot by itself
guarantee faster local Newton convergence.

Keep the existing production default. Preserve the complete experimental
patch and data, but do not add a selectable fallback, molecule-specific
metric rule or empirical weighting parameter. No physical accuracy
threshold was loosened. The default production sources were restored using
only the reverse of this experiment's own patch; the unrelated user changes
were preserved. No default/baseline calculation was rerun after restoration.

`metric.patch` applies cleanly to the parent and includes projector-metric,
bilinear symmetry/positivity, shifted-diagonal, CI-scale and unchanged-cap
tests. Original gradient/HVP derivative checks also passed. `results.csv`
was generated from raw logs by the common `collect.py`.

Local raw records: `results/diagnostics/casscf-metric-20261017/` (ignored).
The remote source hashes, executable hashes, native time reports, input
hashes, node/CPU/affinity, trajectories and qualification logs are retained
under
`/home/guqqgroup/taoxia/work/xmvb-ciah-subspace-20261017/benchmarks/casscf-metric-20261017/`.
