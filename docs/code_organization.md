# Code Organization and Naming

This document defines the production module boundaries and local naming rules.
The directory provides the context, so filenames should describe only the
operation or data owned by that directory.

## Production Modules

```text
src/
  core/     Small numerical and execution utilities shared across methods
  mcscf/    Method-independent CI, optimization, and state contracts
  casscf/   Orthogonal CASCI model and active-integral preparation
  vbscf/    VB structures, nonorthogonal orbitals, derivatives, and workflows
  input/    XMI deck representation and construction of VBSCF inputs
  libcint/  AO integral providers and libcint-specific preparation
  guess/    Initial-orbital and restricted-HF construction
  output/   Molden and accepted-iteration serialization
  cli/      Command-line parsing, execution, and reporting
  tools/    Developer diagnostics and benchmarks
```

There is no `application` or `runtime` layer. The command-line interface is one
consumer of independently named production components, not the owner of their
scientific logic.

## Build Dependencies

The method libraries depend on a common numerical library:

```text
xmvb_vbscf ----> xmvb_mcscf <---- xmvb_casscf
```

An arrow means "depends on". The common library owns direct-CI actions, spin
operations, state weights, matrix-free optimization, and the shared `core/`
eigensolvers. It does not depend on either physical model. CMake rejects direct
VB includes in `mcscf/` or `casscf/`; public solver tests link only
`xmvb_mcscf`, not `xmvb_vbscf`.

The I/O components depend toward the VB model: `xmvb_libcint`,
`xmvb_input_deck`, `xmvb_guess`, `xmvb_input`, and `xmvb_output` form an acyclic
graph. `xmvb_input_deck` is separate because guess construction consumes parsed
deck data, while complete input loading consumes the guess builder.

The `xmvb` executable owns only files in `src/cli/`. Diagnostic executables are
excluded from the default build and live in `src/tools/`.

## Common Solver Boundary

```text
mcscf/
  ci/                          Spin strings and direct-CI actions
  states/weights.*             State-weight validation and normalization
  optimization/
    actions.hpp                Vector action and dual-norm contracts
    neo/                       NEO/CIAH subspaces and numerical solvers
    krylov/                    MINRES and HVP basis operations
    preconditioners/           Algebraic BFGS, coarse, and shifted-metric actions
    trust_region/              Spectral trust-region solve
    trace.hpp                  Scalar solver-event logging
```

These APIs use `xmvb::mcscf`. A `VectorAction` maps a vector to a vector;
the caller defines its coordinates. A `DualNorm` evaluates a covector norm
using the physical model's metric. A common solver must not infer orbital
coordinates, CI normalization, or a VB structure representation from either
callback. `NeoProblem` supplies the gradient and action/metric/preconditioner
callbacks; `CiahWorkspace` owns only numerical subspace state.

VB orbital charts, CI transport, physical curvature construction, and macrostep
evaluation remain in `vbscf/optimization/`. In particular,
`vbscf/optimization/neo/coordinates.*` and `curvature.*` are VB model adapters,
not common solvers. Their dependence on VB derivatives is intentional.

There are no forwarding headers or namespace aliases at the old solver paths.
Callers include the owning module directly. Shared state weights have one
implementation; workflow and derivative code must not normalize independently.

The common integral layer is the next extraction step. AO/active-space kernels
still under `vbscf/integrals/` are not yet method-independent APIs. Likewise,
`casscf/` currently exposes CASCI, not a complete CASSCF orbital optimizer.

## Filenames

- Do not add a `cpp_` prefix; C++ is already the language of the repository.
- Do not encode obsolete implementations with `legacy`, `compat`, or
  `fallback` names. Delete retired paths.
- Avoid umbrella names such as `application`, `runtime`, `manager`, or `utils`
  when a scientific or I/O responsibility can be stated directly.
- Prefer a short noun for a data contract and a short operation name for an
  implementation, using the containing directories to supply context.

## Dimension Names

Use these local abbreviations consistently in numerical kernels:

| Name | Meaning |
| --- | --- |
| `n_bf` | number of AO basis functions |
| `n_ao` | number of active orbitals |
| `n_orb` | total number of orbitals |
| `n_occ` | number of occupied orbitals |
| `n_inact` | number of inactive orbitals |
| `n_det` | number of determinants |
| `n_str` | number of VB structures |

Public data contracts may retain explicit names such as
`n_basis_functions` and `n_active_orbitals`. Short forms are for local
expressions whose scope makes the meaning unambiguous. In particular, never use
`n_ao` for the AO basis dimension; use `n_bf` so it cannot be confused with the
active-orbital count.

Dense physical vectors and matrices are Eigen objects. `std::vector` remains
appropriate for ragged topology, discrete indices, packed integral records,
and raw external I/O buffers, but not as an implicit dense-matrix container.
