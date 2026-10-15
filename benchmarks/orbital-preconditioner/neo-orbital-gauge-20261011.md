# Complete orbital curvature without horizontal structure response

## Experiment

`neo-orbital-gauge.patch` replaces only the core-only NEO branch's orbital
action with `apply_orbital_coupling(vector).orbital_hessian`. This supplies
local-active curvature and overlap-normalization gauge response, without
solving the horizontal structure response equation. Existing admission,
coupled branch, geometry, preconditioning, early trials, trust-radius updates
and final convergence criteria are unchanged. This is an experimental binary,
not a production default change.

Hanhai25 compute job `281651` builds the executable and runs NEGPIA;
array `281652` runs the other five independent systems in parallel. Both
variants of each case run sequentially on the same node, with the same saved
historical input, exact Libcint integrals, Davidson, 32 compact-affinity OpenMP
threads and one BLAS thread. All 12 applications exit zero and explicitly
report convergence. No optimizer iteration cap was added. The projected
gradient and adjacent energy thresholds remain 1e-3 and 1e-7 Eh.

Current binary SHA256:
`71fe190bd423c658e165693b2c925ccc4ccf5cf2498a4054711a947b854d8807`.
Orbital-gauge binary SHA256:
`4cc7c676f06c63821c862762877c019f2c206d742d4bc51b78952ccc578305b1`.
All input hashes match within pairs. Initial energies agree within 1e-12 Eh.
Job `281650` failed patch application before numerical execution and was
replaced by `281651`. The execution mirror source and executable are restored
after building the isolated variant.

## Native results

Each cell is **core-only / complete orbital curvature**. Times are single
matched native samples; small differences are not repeated statistical estimates.
RSS is GNU-time process peak, not requested Slurm memory.

| System | Points | Orbital actions | Rejected trials | Total s | RSS MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| NEGPIA | 17/7 | 168/45 | 4/0 | 20.960/8.822 | 3781.6/3737.3 |
| 240 | 8/8 | 12/13 | 0/0 | 2.182/2.304 | 1772.2/1732.4 |
| 67970 | 7/7 | 13/13 | 0/0 | 3.544/3.646 | 2045.5/2031.5 |
| CONYAH | 7/7 | 13/13 | 0/0 | 1.921/2.020 | 1730.4/1580.6 |
| QOCGUM | 7/7 | 15/14 | 0/0 | 7.555/7.707 | 2036.9/2060.3 |
| YEPKAF | 7/7 | 25/23 | 0/0 | 9.348/9.712 | 2235.4/2234.1 |

Both variants have zero horizontal structure-response actions throughout.
This does not mean zero H/S directional work: the reused coupling API still
constructs overlap/Hamiltonian directional images and coupling forcing before
returning its orbital image. The latter forcing is discarded by this experiment.
No full orbital Hessian or new per-thread integral arrays are introduced.

| System | Core-only final energy, Eh | Complete orbital final energy, Eh |
| --- | ---: | ---: |
| NEGPIA | -436.262505347362 | -436.262507071107 |
| 240 | -343.446302347233 | -343.446302347041 |
| 67970 | -442.247231402783 | -442.247231407475 |
| CONYAH | -343.417822678080 | -343.417822681842 |
| QOCGUM | -480.094766640464 | -480.094766654683 |
| YEPKAF | -1137.242583650810 | -1137.242583617490 |

Endpoint differences are below 3.34e-8 Eh except NEGPIA (1.724e-6 Eh).
NEGPIA's complete-orbital energy is within 1.84e-8 Eh of the preceding coupled
endpoint. Both matched variants satisfy the unchanged stopping criteria;
their endpoint energies are not identical, so no equal-final-energy-accuracy
speedup claim is made.

## Comparison with coupled NEO

The preceding coupled experiment reached NEGPIA in 6 points, 45 orbital HVPs,
293 structure actions and 8.426 s. This complete-orbital variant reaches it in
7 points, 45 orbital actions, zero structure-response actions and 8.822 s.
Those two runs are from different allocations; their approximately 0.4 s
time difference cannot establish a matched performance ranking.

The other five cases previously used 5 to 7 points with coupled NEO, while
this variant uses 7 to 8. It removes their horizontal response solves but does
not recover the coupled iteration reductions. The current core-only method
already needs very few actions on these cases, and additional per-action work
produces modest observed time increases of 2% to 6% in the matched runs.

## Conclusion

For NEGPIA, restoring complete fixed-structure orbital curvature is sufficient
to remove the observed radius-collapse regression: all seven accepted points
retain radius 1 and have no rejected trials. Horizontal structure inversion is
not required to recover this improvement on this trajectory. This is a stronger
result than simply restoring complete coupled NEO, but not a universal speedup.

Before promotion, investigate which directional coupling work is necessary to
obtain the normalization gauge contribution and whether it can be specialized
or reused. The present API intentionally computes both orbital image and
coupling forcing. Removing that overhead needs a mathematically certified
gauge formulation, not another incomplete local-active action. Then repeat
the broader historical panel; do not add molecule-specific exceptions.

Full results and hashes are in `results/history/neo_orbital_gauge_20261011.csv`.
Raw logs are under `results/studies/locality-20261011/neo-admission-281651/`
and `neo-admission-panel-281652/`, mirrored remotely under
`/home/guqqgroup/taoxia/xmvb-runs/native-20261011/`.
