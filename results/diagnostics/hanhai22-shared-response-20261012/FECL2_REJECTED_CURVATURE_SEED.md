# Lowest-curvature startup seed: tested and removed

Hanhai22 short job 1249581 tested one additional transported direction:
the lowest generalized curvature Ritz vector of the previous macro's
cached Q^T H Q and Q^T G Q. The physical orbital and horizontal coefficient
representatives were reprojected into the new chart, with the structure
part aligned to the relaxed root's phase. New Hessian/metric images were
always evaluated at the new point. No old-point operator image or secant
history was used.

The seed replaced one startup residual direction before the existing
three-direction microstep admission. AH capacity, exact operator, physical
metric, micro policies and all convergence tolerances stayed unchanged.
Only one packed orbital direction and one structure coefficient column
were carried between macros; there were no full matrices or thread replicas.

Independent unit tests checked the physical Ritz eigenvector, metric
normalization, no extra operator evaluation during extraction, fresh images
at startup, rejection of a dependent seed, and exact full-space AH roots.
All five unit executables passed in 1249581.

| Same-node anode31 native experiment | Macros | HVPs | Energy / Eh | Wall / s | Peak RSS / MiB |
|---|---:|---:|---:|---:|---:|
| Archived current scalar diagnostic, 1249566 | 13 | 252 | -2181.624464872092 | 9.13 | 504.66 |
| Curvature-seeded startup, 1249581 | 16 | 318 | -2181.624464869361 | 11.15 | 530.24 |

The original input, exact integrals, bound 32-core placement and numerical
policies match. The output passes the existing convergence tests and energy
accuracy, but the experiment worsens actions, macros and native time. It
does not support preserving this startup policy. It also does not disprove
multi-vector recycling with a different, justified selection rule.

The seed type, driver plumbing, workspace option, extraction method and
experimental tests were removed. The previously retained inverse-spectrum
diagnostic and independently exercised phase test remain. No runtime switch
or failed production path is kept. The remote source mirror is restored;
its numerical executable is restored from the verified pre-experiment
1249573 frozen build. Rejected frozen artifacts remain only as evidence.

The next concrete cost target is the observed 214--226 s analytic diagonal
construction in the iron four-step runs, not another threshold adjustment.
Inspect separately the serial core Coulomb projections and the repeated
core/active density graph sweeps before choosing a bounded batch layout.
