# Fe-porphine optimizer benchmark

Purpose: compare conventional L-BFGS (ISCF=5), block-LBFGS (8), and NEO (9)
without changing the orbital initial guess or variational space. This is a
benchmark design, not evidence that NEO is faster.

| Setting | Value |
|---|---|
| Molecule | Neutral Fe-porphine, 37 atoms, 186 electrons |
| States | Singlet, triplet, quintet, separate state-specific calculations |
| Active space | 8 electrons in 11 orbitals; complete spin-adapted VB space |
| Orbitals | OEO; no frozen orbitals or spatial symmetry restriction |
| AO basis | Repository cc-pVDZ, Cartesian in all calculations |
| Integrals | Exact libcint; no RI approximation |
| Structure solver | Davidson |
| Initial orbitals | Shared quintet-ROHF frame preserving its complete occupied space |
| Convergence | Projected gradient infinity norm 1e-3; energy change 1e-7 Eh |
| Resources | Hanhai25 Slurm; independent 32-core jobs, parallel array (up to 45) |
| Safety limits | 1000 iterations; 2 hours per method; incomplete runs explicitly marked |

Geometry and the Fe 3d/4d projector choice originate from the
[public PySCF Fe-porphine example](https://github.com/pyscf/pyscf.github.io/blob/master/examples/mcscf/33-make_init_guess/18-dmet_cas_for_feporph.py).
Ten spherical d projectors are embedded in the complete Cartesian AO space.
Within the ROHF occupation sectors, the two closed and five virtual directions
with the largest projector weights are selected, together with all four singly
occupied orbitals. The remaining 89 closed orbitals form the core. Separate
core, active, and external Fock canonicalization preserves these spaces.
Consequently, the quintet initial CAS contains the ROHF determinant. Its lowest
CI energy must not exceed the ROHF energy. The saved audit checks orthonormality,
core population, and containment of the complete occupied space.

The original DMET-style guess did not preserve that containment: roughly two
reference electrons were outside its core-plus-active space, and the core
included singly occupied weight. Results using that historical frame must not
be interpreted as results from the corrected initial guess. The original
algorithm has been replaced, not retained as an alternate preparation route.

This is **not a literal reproduction** of the spherical-AO, spatial-symmetry
restricted state benchmarks in Kreplin, Knowles, and Werner, *J. Chem. Phys.*
**2019**, *150*, 194106. In particular, the triplet is not constrained to B1g.
Different minima or active-space character must be reported, not counted as
an optimizer speedup. The shared quintet initial frame may have different
initial energies in different spin sectors; only comparisons within a sector
have identical initial energies.

FeCl3 CAS(9e,12o) is deferred until its original geometry is verified. Do not
substitute an assumed Fe–Cl bond length and call it the literature benchmark.

## Run and evaluate

Transfer a clean source snapshot and these files to Hanhai25. From the run
directory, submit `sbatch run.sbatch /absolute/path/to/source-snapshot`.
This job only compiles and prepares the shared initial orbitals. Submit
`sbatch --dependency=afterok:PREP_JOB_ID run-one.sbatch /absolute/path/to/job-PREP_JOB_ID /absolute/path/to/source-snapshot`
to run all nine calculations in parallel; verify current account limits before
submission. Never serialize independent optimizer comparisons by default.
Compilation and orbital preparation run inside their allocation. Source revision
and benchmark script hashes are saved at submission; binary/input hashes,
CPU, affinity, raw output, exit status, GNU-time wall/RSS and Slurm job ID are
retained. A zero exit code alone is not a convergence certificate.

Report, per state and optimizer: convergence status, final energy, final
gradient and energy change, outer steps, orbital/structure actions, wall
seconds, peak RSS in MiB, and state/active-space identity. Retrieve Slurm
accounting with `sacct` after completion. Add a PySCF CASSCF comparison using
the exported `initial_mo.npy`, the same Cartesian basis and spin sector before
making claims relative to an external CASSCF implementation.

To reuse the converged reference without repeating ROHF, run
`prepare.py SOURCE_SNAPSHOT /absolute/path/to/quintet_rohf.chk` inside a fresh
Slurm working directory containing `geometry.xyz`. The checkpoint molecule,
geometry, spin, and AO basis are checked before use. `check-guess.sbatch` prepares
the corrected frame, independently audits it, and evaluates quintet CASCI and
one bounded block-LBFGS step; this is a guess diagnostic, not a convergence run.
