# Direct-CI data layout

The logical levels are VB structures, unique alpha/beta spin strings, and
generated spin-product/channel tiles. An indexed spin product is not an
independent variational coefficient. Davidson vectors remain in the selected
structure space. Storage must follow each contraction, not a single generic
sparse-matrix convention.

## Formats and traversal

| Object | Stored layout | Main consumer |
|---|---|---|
| Original touched product to structure incidence | CSR: product offsets plus one coefficient/structure entry per edge | Expand a Davidson block on original support |
| Reverse structure incidence | Structure offsets plus indices into the same edges | Project images without copying coefficients |
| One-spin single replacements | Fixed-stride, contiguous records | Sigma, integral adjoint, generators and factor moments |
| Transformed support in a beta panel | CSC: beta offsets plus sorted alpha indices | Read/gather compact exact-support columns |
| Multiple RHSs on that support | Dense matrix: support coordinate by RHS | Share a single symbolic pattern across RHSs |
| Hamiltonian channels in a generated batch | Dense row-major channel slab | Contiguous scatter, GEMM, Gram and gather |
| Original-support source factors | Bounded dense one-spin panels | Low-source-count contraction order |
| Exterior maps | Ordered elementary shears and bounded retained row panels | Transform and exact selected projection |

`StructureAction` owns the first two incidence directions; their shared
coefficients must not be copied into a second sparse numerical matrix.
`ColumnPattern`/`ColumnTile` implement the compact CSC panels. They retain
every coordinate in the symbolic support, including zero and tiny values.
The pattern is independent of the current coefficients and is shared across
input/output and RHSs. Full support uses implicit indices, not a stored list
of every row. Cache budgets limit residency, never access to a physical
representation or the selected variational space.

Logical product addresses use `ProductIndex` (signed64-bit), whereas each
single-spin address remains int. Multiplication is widened before packing.
CAS(18,18) already has2363904400 logical products; a selected-space tile
cannot make a signed32-bit global address safe.

## Accepted and tangent ownership

`StructureAdjointState` retains structure coefficients and small integral,
generator and state-gauge tensors. Product order uses the same compact column
format for accepted gradients, orbital directions and their reverse consumers.
No unbounded accepted Y/residual or directional Y/sigma field remains.

Accepted wave-function/residual panels may stay resident under one shared
reservation across all roots, at most one quarter of `--ci-workspace-mib`.
The remaining allowance is divided among input, output, sigma and reverse
workspace; it is not duplicated for each derivative consumer or worker.
A complete resident state avoids replaying its accepted sigma in later HVPs;
coefficient-only responses also borrow its transformed coefficients directly.
Otherwise its primal and tangent images are generated together and immediately
consumed. The representation and numerical accuracy are identical; only data
residency changes. This is not a claim that all accepted work can always be
retained at fixed memory or that streaming has zero recomputation cost.

Exterior rows and their derivatives are shared by input generation and left
projection. A bounded one-body neighbor reader supplies generator cotangents,
so a union of neighboring beta columns cannot silently become a full carrier.
The upper generator cotangent is sufficient for the Cholesky pullback; it is
not a replacement for an oriented transition density. The overlap derivative
uses the exact transpose identity to avoid traversing both operand orders.
The forward orbital direction retains its prepared integral tangent and small
adjoints for the reverse consumer. None of these arrays is copied per worker.

Compact payload caps are not whole-process RSS limits: immutable spin links,
original incidences, active tensors, exterior-row residency and bounded
temporary slabs also coexist. Original touched support can itself be complete.
Complete-space inverse/preconditioner consumers remain a separate migration.

The channel gather's small inverse-index map is a different object: it maps
the current alpha halo and target groups. It is bounded by the batch planner
and repeatedly queried in the hot loop. Replacing it mechanically with a
generic CSR lookup would introduce searches; it is not a complete-product
numerical array. Such a change needs an explicit traversal/reuse analysis.

## Fixed-stride spin links

For a complete fixed-spin sector, every string with `ne` occupied orbitals
has `ne*(norb-ne)` off-diagonal single replacements. The count is independent
of the string, even if only a few VB structures are selected. This follows
from the computational one-spin completion; it is not closure of the selected
VB space. A fixed-stride table therefore needs no CSR row offsets or padding.

`SpinConnections` in `links.hpp` owns a flat `SpinLink` metadata array and a
column-major occupied-orbital Eigen array. A fixed-stride row contains every
link for one string. `singles()` borrows that row through one pointer/count
view; `occupied()` borrows one Eigen column. They allocate nothing and cannot
outlive the shared table. Numerical sigma directions borrow this same table;
equal alpha/beta sectors share it too. Numerical wavefunction and channel
arrays remain Eigen data; the record vector contains graph metadata.

| Field | Type | Meaning |
|---|---|---|
| `source` | `int` | Unchanged spin-string address |
| `pair` | `uint16_t` | Packed symmetric orbital-pair channel |
| `sign` | `double` | Exact +/-1 weight, consumed directly by floating kernels |
| `created` | `uint8_t` | Created orbital in the oriented generator |
| `annihilated` | `uint8_t` | Removed orbital in the oriented generator |
| `occupations` | `uint8_t` | Sorted occupied orbitals |

The current direct-CI contract admits at most63 orbitals. Orbital IDs are at
most62 and the packed pair ID at most2015, so these integer widths are exact.
No CI coefficient, integral or floating-point precision is narrowed. The
created/removed indices are retained explicitly: a symmetric packed pair is
not sufficient for an oriented orbital generator or its adjoint.

Links preserve increasing packed-pair order. Channel slices can use existing
binary searches; no call constructs a temporary link vector or reorders the
logical excitation order. Generation fills disjoint preallocated string rows. It uses a small
stack list of virtual orbitals, not one heap allocation per string.

On the current 32-bit-int platform the link fields occupy16 bytes per single
replacement, plus one byte per occupied orbital. The old nested representation
used24 bytes per replacement, four per occupied orbital, and six vector
headers/allocations per nontrivial string. The retained record table has two
data allocations per distinct nontrivial spin sector, independent of its
string count. This is a storage/locality change, not a reduction of the number of
physical excitation links.

The first contiguous-table prototype used int8 signs. It reduced metadata
further but added integer-to-float conversions in our node-wise kernels and
did not show a runtime benefit. The retained format uses double signs; PySCF
can use integer signs profitably because its vector loops branch once per
link outside the contiguous spectator span. Copying its field widths without
its traversal is not the same optimization.
The record view also replaces the initial prototype's five per-field Eigen
maps with one row view: most hot loops consume address, channel and weight
together. No runtime switch or alternative table format is retained.

## PySCF correspondence

PySCF's `gen_linkstr_index()` also builds fixed-degree string links, and its
C kernels compress addresses, pair IDs and signs into a continuous table.
Its sigma scatters into a bounded dense channel block, calls GEMM, then
spreads the result. It does not form a generic CSR Hamiltonian over all
spin products.

Sources inspected on2026-10-16:

- [Link generation and triangular pair indexing](https://github.com/pyscf/pyscf/blob/master/pyscf/fci/cistring.py)
- [Compressed C link fields](https://github.com/pyscf/pyscf/blob/master/pyscf/lib/mcscf/fci.h)
- [Blocked sigma scatter/GEMM/gather](https://github.com/pyscf/pyscf/blob/master/pyscf/lib/mcscf/fci_contract.c)

PySCF makes beta strings contiguous in its ordinary CI array; our column-major
panels make alpha strings contiguous. Those are transposed storage conventions,
not different physical spaces. Repeated full-array transposes to imitate its
notation would defeat this design. The useful common principle is to keep the
spectator span contiguous and reuse the excitation metadata around GEMM.

## Why not generic BSR everywhere?

BSR helps when occupied blocks are regular and mostly dense. Arbitrary selected
structures generate irregular symbolic support. Fixed rectangular BSR blocks
would retain forbidden zeros again; block size alone does not repair that.
Here, compressed column indices select exact support, while bounded dense
channel/RHS panels supply the block arithmetic. Neither the full Hamiltonian
nor its overlap matrix is stored as CSR/BSR.

## Next ownership migration

Separate three lifetimes: immutable string/incidence topology, accepted-point
integrals/exterior coefficients/state, and one-direction numerical workspace.
The topology should be shared by forward, backward and directional consumers.
At present sigma directions share the spin links, but exterior shear pairs
and topology reconstruction across accepted points still need a separate audit.

The broad-support derivative consumers share compact column patterns,
excitation tables, the generated channel sweep and the physical upper
generator pullback. Output panels retire after their last contributor, not
after an individual temporal response. Fully resident images are accumulated
before the reverse contraction; bounded partial-image consumption remains
exact when the live frontier cannot fit. Input columns have longer lifetimes
because the generator adjoint still reads their neighbors.

The occupation masks and excitation XOR/popcount operations already use
uint64 bit arithmetic. Covalent support packs eight prefix-count comparisons
into one word. Link generation now passes its already-known removed/inserted
orbital labels directly, rather than reconstructing them from two bit scans.
Further bit-level work should avoid repeated conversions and searches, not
introduce a dense `2^norb` address table or replace floating GEMM with bit logic.
