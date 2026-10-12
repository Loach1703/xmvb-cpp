# Third-program porphine triplet calculation

Submitted on Hanhai22, 2026-10-12:

- Build: **1248208**, 32 CPUs, OpenMolcas v26.06.
- Calculation: **1248206**, 32 CPUs, dependency `afterok:1248208`.

The source is a tracked-tree archive from the existing v26.06 checkout on
Hanhai25, commit `8355057f32d65706a35996b5ab07cac2962bb728`. It is compiled
on a Hanhai22 compute node, with GCC 14.1 and MKL 2024.1, no MPI/GA.
Libxc 7 is reused from the existing conda installation; this calculation
does not use a density functional. HDF5 is enabled for AO mapping checks.

Initial build 1248197 failed before configuration because Anaconda's older
libstdc++ shadowed GCC's runtime. The corrected build puts GCC's runtime
first. The calculation's dependency was updated to the corrected job; no
duplicate calculation was submitted.

Build 1248208 subsequently failed because the compute node could not download
libwignernj. The successful Hanhai25 configuration uses bundled v0.6.0, commit
`0c907ee4e726343e333be1d4c8663172da95cbe8`, static, Fortran interface enabled,
baseline Fortran flags (not global Molcas I8 flags). Its tracked source was
transferred offline to Hanhai22. The build script now installs that dependency
locally and supplies its CMake package prefix to OpenMolcas; no download is
needed for this dependency.
Corrected offline build job: **1248254**; calculation **1248206** now depends
on `afterok:1248254`. The dependency was observed compiling successfully.

## Matched physical problem

Fe porphine, neutral, triplet, CAS(8,11), 89 inactive orbitals, 469 Cartesian
cc-pVDZ AOs, no spatial symmetry, no ECP, no RI/Cholesky decomposition, no
relativistic Hamiltonian. The gateway uses checkpoint geometry in bohr
to avoid an additional geometry conversion. The native cc-pVDZ library's
actual primitives/contractions must match the checkpoint before proceeding.

The shared low-energy orbital frame is exported by the same occupied-frame
construction used for PySCF task 1248127_1. Molcas computes its own initial
spin-adapted CI solution; there is no claim of transferring CPP coefficients
into Molcas's differently ordered CSF representation.

Conversion checks atom matching, Cartesian AO permutation and contraction
phase, primitive exponents/coefficient ratios, normalized overlap, complete
MO orthonormality, and fixed-width INPORB serialization. Any mismatch aborts
before orbital optimization; it is not replaced by a new automatic guess.
Gateway's GUESSORB is used only for basis/overlap metadata, not the final
initial orbitals supplied to RASSCF.

## Optimizer and reporting

Molcas uses RASSCF CASSCF with its super-CI/quasi-Newton orbital optimizer
(`QUNE`), not an asserted implementation-identical copy of PySCF AH or CPP
NEO. The active variational space and starting occupied frame are matched.
RASSCF uses `Thrs = 1e-7 1e-3 1e-3`, `Tight = 1e-10 1e-4`, and iteration
limits 1000/200. Program-specific residual norms are not presumed identical.
There is no short diagnostic time or macrostep cap; Slurm's 15-day partition
maximum applies.

Setup and orbital optimization have separate `/usr/bin/time -v` records.
Use both when reporting full workflow time; report optimization time separately
if useful. Output, conversion certificate, RasOrb/checkpoint files, hashes,
and scratch path are retained under:

`/home/guqqgroup/xiatao/xmvb-runs/shared-response-20261012/porphine-molcas-1248206`.

Read the actual initial/final energies, normal convergence flag, orbital
iterations, elapsed time, and peak process RSS after the job finishes. This
submission is not yet evidence of agreement or performance.
