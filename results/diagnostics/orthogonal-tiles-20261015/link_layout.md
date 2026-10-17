# Direct-CI layout audit and contiguous one-spin links

Parent: `5e8205c4`. Builds and all numerical runs use Slurm compute nodes on
Hanhai22/anode26, eight bound OpenMP workers, BLAS one, unchanged Release flags
and physical tolerances. Historical binaries were not rerun. Other workloads
share the node. Full source/binary hashes and raw native resource logs are
retained. This round is a data-layout migration, not an optimizer change or
a claim of lower Hamiltonian arithmetic scaling.

## Chosen storage contract

The developer contract is `src/mcscf/ci/README.md`. Structure incidence uses
CSR and reverse indices into the same coefficient data. Transformed support
uses compact CSC by beta column. Hamiltonian channels remain bounded dense
GEMM slabs. Complete one-spin sectors have fixed-degree links, so the new
`SpinConnections` now stores a fixed-stride record array and one occupied-
orbital Eigen array rather than nested per-string vectors. A borrowed row
view replaces five per-field maps. Records preserve the former packed-pair
ordering; coefficient precision, structure space and Slater--Condon graph
are unchanged. Sigma directions and equal alpha/beta sectors share the table.

The integer fields are exact within the existing 63-orbital contract: pair IDs
fit uint16 and orbital IDs fit uint8. Signs retain double, the numerical
kernel's direct input type. Source string IDs retain their existing int type.
No tolerance or amplitude pruning is introduced. Each distinct nontrivial spin
table retains two data arrays, independent of its string count. Occupied and
link rows are borrowed views, not newly allocated vectors.

PySCF primary sources inspected on2026-10-16 are linked in the developer note.
They show fixed-degree compressed links and blocked scatter/GEMM/gather,
not a generic CSR Hamiltonian over the spin-product space. Our record layout
retains both oriented orbital indices for generators as well as packed pairs.

## Exact metadata payload and discarded prototypes

On this platform the old retained payload, computed from its source types and
explicit reserved capacities, is U*(24*ne*(norb-ne)+4*ne+144) bytes. It includes
the six per-string vector headers, not the fixed-size owning object. This
old value is a source-derived calculation, not an allocation-profile reading.
The first SoA prototype used U*(9*ne*(norb-ne)+ne) bytes with integer signs.
The retained record layout uses U*(16*ne*(norb-ne)+ne) bytes on this platform.
Temporary construction maps and whole-process RSS are different quantities.

| One spin sector | Strings | Old calculated bytes | Integer-SoA measured bytes |
|---|---:|---:|---:|
| 12 orbitals, 6 electrons | 924 | 953568 | 304920 |
| 14 orbitals, 7 electrons | 3432 | 4626336 | 1537536 |

The integer-SoA 14-orbital table is 66.8% smaller. Its one captured construction takes
0.002142 s, versus saved0.006847 s in1261729. Do not extrapolate this small
initialization saving to total SCF time.

The integer-sign SoA did not show a runtime benefit. It also introduces
integer-to-float conversions in node-wise consumers, unlike PySCF's
sign-branch-outside-the-spectator-loop traversal. We tested double signs in
SoA (1261824/1261835), then retained a simpler AoS record row (1261842):
all consumers read one record view instead of five field maps. These are
exact data-layout alternatives, not numerical thresholds. Neither SoA
implementation is retained as a production option.

The retained record table measures2714712 bytes for the 3432-string sector,
41.3% less than the original nested vectors, with two data allocations instead
of O(U) allocations. Job1261842 passes11 tests; its independently checked
response takes3.67168 s, accepted preparation3.64118 s, process RSS411.965 MiB,
image error1.11e-15, adjoint error6.45e-16. The integer-SoA and double-SoA
response captures were3.82875 and3.94234 s. This single-run series does not
prove causal attribution to a field type or a total-time improvement over
saved1261729 at3.33289 s. Retained metadata is smaller and ownership simpler;
a universal runtime win is not established by this round.

Finally, generation passes already-known removed/inserted labels directly
instead of decoding them from the XOR mask twice. It still forms the source
mask by XOR and the fermionic sign by popcount. This removes unnecessary bit
work, not the physical excitation graph. Job1261850 compiles and validates
this final change together with the retained record layout.

## Integer-SoA prototype: validation and response timings

Jobs1261793 and1261799 each pass the11 focused tests. Additional tests verify
the exact shared-table payload and empty/full-occupied spin sectors. These
exercise zero-link borrowed maps without inventing a separate kernel path.
The shared changes cover sigma, sigma/integral adjoints, spin factors, oriented
generators, projected structure responses and their existing derivative tests.

Same CAS14/14 selected132 fixture, 3432 strings per spin:

| Quantity | Saved1261729 | New1261799 |
|---|---:|---:|
| Accepted preparation / s | 3.17182 | 3.01297 |
| Cold response / s | 3.33289 | 3.82875 |
| Whole native test process / s | 7.56 | 7.76 |
| Native process peak RSS / MiB | 438.145 | 389.832 |
| First sweep prepare / s | 0.42203 | 0.40792 |
| First sweep read / s | 0.95669 | 1.16169 |
| First sweep compute / s | 1.43049 | 1.71467 |
| First sweep consume / s | 0.32939 | 0.35389 |
| Input/output cache evictions | 0/0 | 0/0 |

The new numerical checksum matches the saved printed value. Independent
complete-carrier image and adjoint errors are1.11e-15 and6.45e-16. Generated
channel counts and numerical slab sizes are unchanged. The compute timer
includes scatter, GEMM and Gram; it is not an isolated BLAS measurement.
Affinity differs: saved CPUs152-159, new136-143. The measurements do not
establish a total-time speedup or attribute the regression to either sign
conversion or external scheduling. In particular, the process-RSS change
cannot be attributed entirely to the much smaller metadata difference.

## Integer-SoA prototype: complete SCFs

Job1261805 passes11 tests and completes nine240 optimizations: selected1/3
with dense/Davidson and COV/ION/FULL. All step/HVP counts are retained, and
final energies agree with the saved baseline within3e-12 Eh. Compared with
the immediate saved1261730 family SCFs, this capture is slower:

| Space | Steps/HVP | Saved SCF / s | New SCF / s | New RSS / MiB |
|---|---:|---:|---:|---:|
| COV | 3/20 | 6.3685 | 6.5086 | 1529.44 |
| ION(0:1) | 3/27 | 8.5650 | 9.1902 | 1515.07 |
| ION(0:2) | 4/32 | 10.0897 | 11.2736 | 1520.41 |
| ION(0:3) | 4/32 | 10.5578 | 11.3475 | 1510.83 |
| FULL | 4/30 | 8.9645 | 10.2447 | 1518.66 |

Job1261807 passes the eight additional geometry/optimizer tests and converges
FeCl2:15 steps,247 HVP, -2181.624464846785 Eh,11.0093 s SCF,410.34 MiB RSS.
The immediate saved1261731 was16 steps/276 HVP, -2181.624464870515 Eh,
13.2075 s,415.74 MiB. The final energy difference is2.37e-8 Eh under unchanged
outer criteria. These changed action counts are not evidence of a new
convergence method or per-action speedup; do not explain them as rounding
without further analysis.

Raw parsed tables: `link_layout_responses.csv`, `link_layout_240.csv`,
`link_layout_selected_240.csv`, and `link_layout_fecl2.csv`.

## Retained record layout: final SCF qualification

Job1261850 passes11 tests and all nine240 optimizations after removing the
redundant bit decoding. The tests also retain the exact shared-layout byte
checks for unequal, empty and full-occupied spin sectors. No optimizer strategy
or numerical threshold changes are included.

| Space | Steps/HVP | Final SCF / s | Final RSS / MiB | Final energy / Eh |
|---|---:|---:|---:|---:|
| COV | 3/20 | 6.6855 | 1520.36 | -343.239211501260 |
| ION(0:1) | 3/27 | 8.6599 | 1510.54 | -343.407736636065 |
| ION(0:2) | 4/32 | 10.7693 | 1547.23 | -343.443415326851 |
| ION(0:3) | 4/32 | 11.7334 | 1515.35 | -343.446256564845 |
| FULL | 4/30 | 9.7405 | 1515.57 | -343.446302545200 |

Step/HVP counts match the saved baseline; energies differ by at most2e-12 Eh.
Selected1/3 dense/Davidson controls all take3 steps/21 HVP and agree in energy
within1e-12 Eh. The immediate pre-layout1261730 was faster on most family
captures; this round is not a universal end-to-end acceleration.

Job1261856 passes the eight additional tests and converges FeCl2 in16 steps,
281 HVP,12.4227 s SCF,402.83 MiB RSS, final energy -2181.624464871859 Eh.
The energy differs from saved1261731 by1.34e-9 Eh. Although wall time is lower
than its13.2075 s, the HVP count increases from276 to281; no convergence
improvement is inferred from a representation-only change.

Final tables are `link_layout_final_240.csv`,
`link_layout_final_selected_240.csv`, and `link_layout_final_fecl2.csv`.

## Bit-operation audit

Occupation masks, XOR-generated replacements and popcount signs were already
present. `CovalentSupport` performs eight byte-lane prefix comparisons in one
word, under its documented no-carry bounds. The removed/inserted orbital
labels are now passed directly from the enumeration loop, eliminating two
bit decodings that reconstructed already-known information.

Static linked-binary inspection runs on a Slurm compute node (job1261876);
no compiler flags or workload are changed. The retained assembly shows the
sign parity calculation folded by shifts/XOR into an x86 `cmovp`, at addresses
0x8fe11c--0x8fe139. It does not compute an unnecessary full population count.
This establishes the instruction sequence, not its wall-time share. The first
object-level scan1261869 found no matching mnemonics and exited1. The next
scan1261873 selected mangled symbols while requesting demangled matching and
produced empty function bodies; it is not valid ISA evidence. The corrected
scan1261876 retains1026 lines of actual linked code. Its initial mnemonic
filter did not include `cmovp`; the full assembly supplies that evidence and
the checked-in filter now includes parity instructions. No inference of
missing optimization is made from a missing POPCNT mnemonic alone.

## Next concrete issues

The compact layout removes scattered metadata ownership and provides the
shared format required by directional tiles; a universal runtime improvement
is not demonstrated by this round. Further kernel work must exploit the
contiguous spectator span, rather than assume narrower integer fields suffice.
Broad-support accepted/directional floating-point product arrays remain.

An independent correctness issue found in this audit is signed32-bit flattened
product IDs. CAS18/18 has48620 strings/spin, whose2363904400 logical products
exceed the range even for a tiny touched support. The root TODO records the
exact source locations. Fix the product-ID contract and add an index-only
overflow test before claiming that larger selected spaces are supported.
