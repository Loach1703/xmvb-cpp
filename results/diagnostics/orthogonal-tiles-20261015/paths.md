# Selected-structure prefix-path verification

Date: 2026-10-15. Production H/S scheduling and derivative storage are not
switched to this prototype. These are path-generation timings, not molecular
VBSCF timings or end-to-end speedups.

## Provenance

- Hanhai22, Slurm `1259298`, `anode62`, AMD EPYC 9654.
- Eight allocated CPUs, explicit `srun --cpus-per-task=8 --cpu-bind=cores`;
  affinity `118-125`. **The path generator is serial.**
- GCC 14.1.0, existing Release build and dependency environment; exact binary,
  tested source hashes, CMake cache, CPU details and native GNU-time logs are
  in `paths-1259298/`.
- Final extra leave-and-return assertion: Slurm `1259316`; logs and hashes in
  `paths-validation-1259316/`. Both jobs passed 12/12 relevant CTests.
- Base source revision `07f7e883` plus the selected-support changes and the
  new test prototype. No historical executable was rerun.

## CAS(12,12) support generation

The selected structures are the first N noncrossing Rumer pairings in the
test's fixed order. Different subsets of the same size can have different
support sizes. The full product universe has 924 x 924 = 853,776 entries.

| Structures | Original touched products K | S products | I products | Build (s) | Generate S (s) | Generate I (s) | Metadata + traversal (KiB) | Process peak RSS (KiB) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 64 | 136374 | 523307 | 0.0000314 | 0.005347 | 0.027141 | 2.742 | 6308 |
| 3 | 128 | 157108 | 559367 | 0.0000506 | 0.006291 | 0.034528 | 5.484 | 6332 |
| 20 | 376 | 192380 | 603701 | 0.0000931 | 0.010585 | 0.063692 | 16.453 | 6408 |
| 132 | 924 | 226512 | 608751 | 0.0001210 | 0.017427 | 0.164182 | 41.133 | 6544 |

All emitted counts equal the previously saved exhaustive support counts.
The generator allocates no product-space mask or carrier. Its metadata and
traversal sizes above exclude seed-input construction and the test fixture;
GNU-time peak RSS includes the whole test process. An I traversal also marks
S membership, so production need not generate both separately just to label
paths. The 132-structure test deliberately uses the general bitset rule, not
the full-covalent closed-form shortcut.

The generic bookkeeping is measurable: I grows by only 16% from the first
to the last row, while the number of seed words grows from 1 to 15. Thus
direct enumeration is not evidence that all support overhead is negligible.
The algebraic full-covalent simplification removes that seed-word factor.

## Correctness

Small-space tests compare every generated product against an independent
sorted-occupation dominance oracle and explicit one-electron neighbors.
They reject duplicates and verify the S label on every emitted I product.
Cases include ionic and mixed seeds, unequal spins, mismatched-seed traps,
duplicate/empty seeds, a zero-electron spin sector and bit index 62.

A 30-orbital, 15+15-electron **ionic product** with both spins in the lowest
15 orbitals generates S=1 and I=451 without constructing either complete
single-spin space. This is a large-universe enumeration test, not a
30-orbital covalent-VB molecular calculation.

The ordered-generator Hamiltonian contraction sums only generated I paths
and is compared with the packed direct-CI result:

| Orbitals, alpha, beta | Maximum absolute bilinear error | Contribution omitted if I is truncated to S |
| --- | ---: | ---: |
| 4, 2, 2 | 3.55e-14 | 0.0169055 |
| 5, 3, 2 | 2.27e-13 | 0.0399336 |
| 6, 3, 3 | 1.48e-12 | 0.1470253 |
| 5, 2, 1 | 2.84e-14 | 0.0258194 |

These are unnormalized synthetic bilinear forms, not molecular energies.
The two contractions share one-body generator links; the comparison tests
support/ordered-versus-packed contraction consistency, not an independent
fermionic-sign implementation. Dense generator images exist only in this
small oracle. They are not the proposed production workspace.

## Next implementation boundary

The next step is a tile scheduler that consumes the generated intermediate
products and their single-generator links directly. Existing ambient-beta
loops, transform-panel regeneration and derivative carrier allocation are
not removed by this generator alone. Measure those consumers after migration;
do not infer a total-speed benefit from the support table above.
