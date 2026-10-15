# Orthogonal-CI structure actions

## Scope

`StructureAction::apply`, `apply_overlap`, and `apply_streamed` use the same
panel pipeline whenever both spin-string sets are complete. Davidson vectors
remain in the selected VB structure space. Completeness is a representation
condition, not a memory-admission test. No VB structures are added by the
change of representation.

The dense-carrier APIs used by selected-state adjoints and orbital derivatives
are separate consumers. This change does not claim to remove their allocations.

## Data flow

1. Generate a small set of rows of the beta exterior transform by applying
   its adjoint to coordinate rows. Never build the full compound matrix.
2. Contract the sparse structure expansion directly with these rows, then
   apply the alpha exterior transform. This produces one orthogonal column
   panel without a full raw product array.
3. The shared direct-CI sigma requests beta columns through a reader. Its
   occupancy and single-excitation links scatter both spins into a shared
   one-body intermediate. The existing absorbed Hamiltonian pair kernel
   performs the same contraction as the materialized sigma implementation.
4. Sigma contributes columns through an additive writer. Contributions are
   accumulated in resident output panels, not complete Hamiltonian images.
5. On consumption, apply the alpha exterior adjoint and sample the beta
   transform only at the sparse expansion's touched spin products. Accumulate
   directly into structure-space output. A touched product's dot product is
   evaluated once and shared by all structures containing it.

Input panels have an LRU residency bound. An evicted input is regenerated
exactly if requested again. An output panel is projected before eviction;
later contributions begin a new zero panel and are projected additively.
Thus eviction changes execution order, not the operator. Overlap contributions
are consumed once per input panel, even when that panel is regenerated.

## Storage and reuse

- Structure RHS batching is bounded by the sigma scratch size. Beta-transform
  rows are shared across RHS, input expansion, and output projection.
- The default input and output panel caches each target 64 MiB. At least
  one column must fit as an unavoidable linear-in-string workspace.
- Transform rows have their own bounded panel cache, also targeting 64 MiB;
  entries are panel-width times the beta string count, never an unbounded
  compound matrix.
- Sigma gathers only the alpha single-excitation halo of its current target
  tile. Beta columns are distributed over workers, which reuse their small
  scatter/GEMM/reduction tiles. Shared halos and all worker buffers are sized
  together under a 64 MiB target; there is no carrier-sized copy per thread.
  Buffers are reused across alpha sweeps instead of reallocated for each tile.
- One persistent OpenMP team covers all sigma sweeps, following PySCF's
  `FCIcontract_2e_spin1` organization. The active tile's excitation halo is
  gathered once, beta columns are distributed over workers, and buffered
  contributions are reduced before the next tile. Callbacks execute serially;
  callback exceptions are rethrown after the team exits. Resident source panels
  are prepared before this team starts, retaining transform parallelism.
- There is no unbounded product-grid allocation in this H/S pipeline and no
  unbounded compound-transform allocation. Small carriers may naturally have
  all their panels resident; large ones are evicted under the same policy.
- `apply_streamed` accepts a workspace target to control panel residency and
  pair-tile size. This target is not a bound on all solver/integral memory and
  never changes the selected structure space or rejects a complete carrier.

The budget controls reuse, not a proof of computational tractability. A very
large carrier can still incur substantial regeneration and sigma work.

## Validation

`test_structure_tiles` compares the pipeline with the independently assembled
materialized reference. Cases include complete covalent Rumer spaces,
selected structures, unequal alpha/beta sectors, multiple RHS, and forced
one-column eviction. The shared sigma callback also has a multiple-RHS parity
test.

The default test also includes a 252-by-210 unequal-spin carrier with two
RHS, covering multiple alpha sweeps and a short final tile, and verifies
reader/writer exception propagation.

`test_structure_tiles --benchmark12` uses the 132 independent singlet covalent
Rumer structures for twelve singly occupied orbitals: 924 strings per spin,
924 touched raw products, and 8448 structure-expansion terms. It is a synthetic
fixed-orbital algebra test, not a molecular SCF benchmark.

`--benchmark12-no-reference` excludes the reference carrier from RSS
measurement. `tests/mcscf/tiles.sbatch` additionally runs the actual LOFLEA
CAS(12,12) input with exact integrals, NEO and structure-space Davidson.

### Persistent-team validation (Hanhai22 job 1257853)

All five targeted CTests passed. The 32-core synthetic test used 924 strings
per spin and 132 covalent structures. Hamiltonian and overlap absolute Frobenius errors
against the materialized reference were 1.24e-14 and 2.63e-15. Native process
peak RSS was 92,120 KiB, excluding the materialized reference. First and warm
actions took 0.227 and 0.125 s; the parity process measured a warm tile action
of 0.111 s against a warm materialized reference of 0.065 s. This verifies
equivalence and bounded execution, not a speedup. In particular, the current
streaming overhead remains significant for a carrier that fits in memory.

Reference execution organization:
[PySCF fci_contract.c](https://github.com/pyscf/pyscf/blob/master/pyscf/lib/mcscf/fci_contract.c).
PySCF itself retains complete coefficient/image arrays; the sparse structure
expansion, exterior-transform panels and exact selected-space projection here
are additional VB-specific operations.
