# Coupled inverse diagnostic and iron four-step accounting

## Cached-image experiment: Hanhai22 1249573

All five unit executables passed, including the independent positive/signed
inverse reference and a check that diagnostic spectra evaluate no HVPs.
The formerly defined but uncalled CI restart-phase test is now called by
`test_ciah`; its actual execution passed too.

For cached directions Q and shifted action L=H-omega G, inspect
eigenvalues of (LQ)^T P(LQ) relative to Q^T LQ only when Q^T LQ is positive
definite. The physical curvature spectrum uses Q^T GQ. Only matrices of at
most 30-by-30 and one inverse image are allocated. No complete Hessian,
metric, additional HVP or per-thread carrier is constructed. The skew of
the inverse bilinear form is reported before symmetrization.

| Macro | Raw physical curvature range | Base inverse spectrum | Balanced inverse spectrum |
|---|---|---|---|
| 9 | 0.00560--170.56 | 0.0818--20.71 | 1.0000--1.0000 |
| 10 | 0.00528--193.35 | 0.0369--37.60 | 0.9993--1.0027 |
| 11 | 0.00524--826.81 | 0.0520--43.91 | 0.9993--1.0148 |
| 12 | 0.01238--226.29 | 0.0333--45.78 | 0.9998--1.0679 |

Inverse bilinear-form relative skews are about 2e-16--5e-16. First two
macros have indefinite projected curvature; no SPD inverse condition
estimate is assigned to them. Native event ordering, all policies and
acceptance remain unchanged.

These are trial-space Ritz measurements, not the spectrum of the entire
operator. The near-unit balanced spectrum is expected from exact coarse
balancing on the covered sample space. It does not certify a good inverse
on the complement. The new evidence rules out blaming an inaccurate coarse
inverse on already represented modes, not all preconditioning issues.
Together with the earlier keyframe projection data, the unresolved work is
learning missing modes and not relearning useful modes at every macro reset.

The audited run converged in 12 macros/247 HVPs at -2181.624464872044 Eh,
10.49 s and 490.21 MiB GNU RSS. Full diagnostic mode includes additional
connection and inverse evaluations, so its wall time is not a speedup
measurement. It must not be compared as a native optimization variant.

## Iron porphine: completed bounded four-step experiment 1249554

Both new tasks ran on anode07, 32 bound cores, exact integrals, original
input orbitals; application exit 2 means the intentional four-macro limit,
not convergence. Archived 1249418 runs are not repeated. Their nodes differ
(singlet anode31, triplet anode45), so wall differences are exploratory.

| Spin | Archived/new HVPs | Archived/new fourth energy / Eh | Archived/new wall / s | New GNU peak RSS / GiB |
|---|---|---|---|---|
| Singlet | 36/51 | -2244.717977952240 / -2244.719899599181 | 424.92/433.25 | 38.622 |
| Triplet | 42/55 | -2244.769462512970 / -2244.759927208631 | 621.89/619.33 | 38.632 |

There is no demonstrated iron speedup or universal trajectory improvement.
In particular, the new triplet fourth-step energy is higher. Do not start
a blind long run from this evidence or claim comparison to the converged
PySCF energies has succeeded.

| Timed section / s | Singlet | Triplet |
|---|---:|---:|
| Analytic orbital diagonal construction | 214.20 | 225.59 |
| Expansion HVPs, including image observer | 27.14 | 29.44 |
| Complete structure norms | 2.79 | 2.78 |
| Forward CI solves | 10.95 | 167.18 |
| Forward active two-electron transformation | 15.81 | 22.71 |

The analytic diagonal builds perform separate density-operator sweeps for
each core orbital and each exact active-density-span column. With 89 core
orbitals and 23 active source columns, this is 112 density sweeps per build,
in addition to pair products and frame contractions. This is not 112 HVPs
and must not be hidden in their counter. A bounded shared batch contraction
is a separate performance opportunity; its numerical equivalence and memory
budget need validation before changing production. Changing the diagonal
precision or silently replacing it with the old one-electron model is not
an acceptable substitute.
