# Porphine triplet: matched low-branch convergence

## Submission

Hanhai22 Slurm array **1248127**, submitted on 2026-10-12:

| Task | Method | CPUs | Input |
|---|---|---:|---|
| 1248127_0 | CPP NEO/CIAH, analytic-diagonal preconditioner | 32 | Exported triplet endpoint from 1247942 |
| 1248127_1 | PySCF default one-step CASSCF AH | 32 | Same occupied frame and actual spin-string carrier from 1247963 |

Both start from the same physical wavefunction used in the completed short
cross-program restart, not the old high-energy PySCF endpoint or either
two-macrostep endpoint. The CPP input retains CAS(8,11), 89 inactive orbitals,
469 Cartesian cc-pVDZ AOs, OEO, full triplet structure space, exact libcint
integrals, no ECP, and the existing 1e-7 Eh / 1e-3 gradient thresholds.
PySCF retains its previously tested settings and gradient convention.

The CPP binary is frozen from final-source regression 1248081. Input, basis,
and executable hashes are saved in the run directory. CPP does not rebuild
or change optimizer numerics for this run. PySCF skips the extra independent
directional-gradient/AO2MO diagnostic used in the earlier short audit.

Both use 1000 macro iterations as the existing safety limit, not a short
diagnostic cap. Slurm's partition maximum is 15 days; this is a scheduler
limit, not an algorithmic stopping criterion. There is no custom short wall
limit. No extra convergence or endpoint stability criterion is introduced.

PySCF writes `cross-restart.chk` during optimization. CPP retains the input
endpoint, scalar trace, full output, and final Molden output on normal exit;
it does not currently provide an equivalent periodic orbital checkpoint.

Remote artifacts:
`/home/guqqgroup/xiatao/xmvb-runs/shared-response-20261012/porphine-low-branch-1248127/{neo,pyscf-ah}`.

At initial inspection both tasks were pending due to Priority. Results must
be reported only after reading the actual outputs and Slurm accounting.
Compare initial/final energies, normal convergence, macro counts, CPP HVP
counts, available PySCF action counts, wall time, and peak process RSS.
Full process timings include each program's integral/setup/output work;
do not relabel them optimizer-only timings.

## Remaining electronic offset: source inspection

The existing transferred-start offset is approximately 1.681972e-6 Eh.
Basis/ECP/overlap and nuclear repulsion have already been checked. No new
fixed-point energy comparison was submitted.

The loader calls `LibcintMaterializedIntegralProvider::build` with default
options. These enable historical shell-pair prescreening and discard
normalized individual AO ERIs below 1e-10. The prescreen uses the geometric
mean of raw shell-pair maxima. Consequently the materialized AO graph does
not contain every libcint ERI, even though the path is called exact rather
than RI. This is a candidate for the accumulated electronic offset, not a
demonstrated explanation. AO normalization and contraction multiplicities
also remain to be distinguished by integral-level checks.

No production screening, integral tolerance, or optimizer threshold was
changed. Next diagnose matrix/block differences at the integral level,
without another redundant fixed-wavefunction energy comparison or a
complete four-index AO allocation.
