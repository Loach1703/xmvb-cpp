# Orthogonal CASSCF

Use `CASSCF` instead of `VBSCF` in `$CTRL` of the existing `.xmi` format.
`NAO`, `NAE`, `NMUL`, `BASIS`, `INT=LIBCINT`, `GUESS`, `$GUS`, `$GEO` and
`ITMAX` retain their input meanings. `NAO` is the number of active orbitals,
not AO basis functions. The molecular charge sets the total electron count;
the remaining electrons form doubly occupied inactive orbitals.

The initial implementation is exact-integral, real, single-state, complete
CAS, with pure-spin Davidson CI and the common CIAH optimizer. `STR`, `$STR`
and orbital localization do not restrict the CASSCF variational space.
The supplied core span and active span projected off the core are preserved
when forming a full AO-metric orthonormal frame. Thereafter updates are finite
antisymmetric MO rotations, not the VB nonorthogonal coefficient update.
`ISCF` continues to select VB optimizers for VBSCF; this CASSCF implementation
always uses CIAH (`ISCF=9` is the natural explicit choice).

RI, state averaging, dense CI and other orbital optimizers are not enabled in
this first implementation. Unsupported explicit selections are rejected.
Converged `.orb` files use the existing fixed-point `$GUS` body format without
block delimiters. `MOLDEN` uses the same basis-and-orbital output writer as VB.

## Ownership

- `orbitals`: nonredundant rotation coordinates and one-time frame completion.
- `integrals`: one exact accepted-point AO→CAS preparation.
- `derivatives`: fixed-state energy/gradients and fused coupled actions.
- `driver`: CI relaxation and physical macro acceptance.
- `mcscf/ci`, `mcscf/optimization`, `integrals`, `linalg`: shared kernels.
- `input/loading/casscf.*`, `cli/casscf.*`: `.xmi` and output adapters.

No physical CASSCF source includes a VB objective, chart or structure type.
The input/output adapters reuse the existing guess and orbital-file contracts;
they do not generate VB structures. No orbital Hessian or CI Hamiltonian is
assembled. Complete-CAS CI vectors still scale as the alpha/beta string
product; bounded sigma tiles do not remove that intrinsic CI storage.
The analytic diagonal currently requires the shared AO CSR representation.

The coupled preconditioner reuses `mcscf::CoarsePreconditioner`: orbital and
CI samples come from the same physical action, with no extra HVP. Its images
are local to one accepted point. Across macros only the last accepted search
direction is recycled, with the CI sign aligned and its tangent reprojected;
old Hessian images are never carried to a different orbital point.

Derivation: `article/orthogonal_casscf.md`. Same-initial-frame Slurm/PySCF
checks: `benchmarks/casscf-20261017/`.
