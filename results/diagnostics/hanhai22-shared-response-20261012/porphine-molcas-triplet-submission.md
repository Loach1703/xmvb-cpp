# Iron porphine: matched Molcas triplet submission

Submitted on Hanhai22 on 2026-10-14 as Slurm job **1253762**, with 32 CPUs
and 192 GiB requested memory. The singlet job **1253591** is unchanged.

The preparation and submission scripts now accept an explicit `singlet` or
`triplet` argument. Both states use the same reference contraction columns,
written to the custom `XMVBREF` basis. This avoids the primitive/contraction
mismatch that stopped the previous triplet conversion before optimization.

- Molecule: neutral FeC20H12N4, 186 electrons.
- Basis: Cartesian cc-pVDZ, 469 AO functions; all-electron, nonrelativistic.
- Active space: CAS(8,11), with 89 inactive orbitals, C1 symmetry.
- RASSCF multiplicity: `Spin = 3`; no RI or Cholesky decomposition (`NOCD`).
- Initial orbitals: original `inputs/triplet/input.xmi` `$GUS`, exported
  without orbital optimization using `export-frame converge original triplet`.
- The converted orbital occupations have 89 inactive double occupations,
  three active double occupations and two active single occupations.
- Conversion checks: normalized AO overlap and full orbital Gram matrix,
  each required to agree within `1e-8`.
- Optimization and CI thresholds are identical to the submitted singlet.

Submission command (from the remote run root):

```sh
sbatch --parsable --job-name=porphine-molcas-triplet \
  porphine-molcas-matched.sbatch triplet
```

Run directory:
`/home/guqqgroup/xiatao/xmvb-runs/shared-response-20261012/porphine-molcas-triplet-1253762`.
Preparation, overlap-only setup, conversion, and optimization have separate
logs. Conventional four-center integral generation precedes RASSCF, so a
running job is not itself evidence that orbital optimization has started.

## Launch validation

The job started on `anode01`. Preparation, overlap-only setup, original-guess
export, and conversion all completed successfully. Measured conversion errors:

- Normalized AO overlap: `1.3655743202889425e-13`.
- Full orbital Gram matrix: `1.6260892464593598e-12`.

The generated RASSCF deck contains `Spin = 3`. The main optimization step is
running in SEWARD conventional integral generation; no RASSCF energy is
available yet. Conversion/export/optimization error logs are empty at launch.
