# Joint orbital-direction sigma contraction

## Scope and algebra

The implementation combines `H Z + deltaH Y` and the integral adjoint in
one generated-target sweep. It removes a second graph/halo traversal, a
second accepted-state scatter, a separate output gather, and a separate
generator-adjoint contraction. The two different channel-kernel GEMMs are
still necessary. Three shared channel slabs replace two, without increasing
the existing 64 MiB numerical workspace budget. This is a constant-factor
optimization, not a new asymptotic scaling claim.

The selected structure space need not be closed under orbital transforms.
The projection, accepted left-transform term, state weights, and certified
spin-exchange multiplicities are retained. No optimizer threshold changes.
The derivation is in `article/source_column_hamiltonian.tex`, Section 4.3.

## Validation

- Hanhai22 job 1262578: all 11 focused CTests passed. New direct-CI tests
  include unequal alpha/beta spaces, widths 1 and 2, exchange parity 0/+1/-1,
  nonclosed support, zero inputs/directions, and optional forcing diagnostics.
- Job 1262587: CAS(12,12), 1/3/20/132 selected covalent structures, current
  action and structure response against the complete-carrier oracle. Maximum
  Hamiltonian error 1.95e-14, overlap error 2.02e-15, integral-adjoint error
  2.78e-15 (reported norm errors).
- Job 1262574: orbital HVP directional finite-difference errors at most
  1.57e-10; energy directional-derivative error at most 6.50e-11. These tests
  use the production low-memory sweep; finite differences are validation only.

## Benchmark definition

The primary baseline is the unmodified GitHub production source at
`946934570c096ba5072923effa7fb063736731f1`. The user authorized one matched
baseline campaign because no saved CAS12 timing exists for that revision.
Job 1262577 failed while compiling the external benchmark driver (a header
moved after the baseline revision); no baseline workload ran. Job 1262586
reuses that library build and fixes only the external fixture driver.

Both variants use the same constructed 12-orbital/12-electron integral fixture,
1/3/20/132 leading Rumer covalent structures, and all 924 strings per spin.
The complete spin tables are provided by the fixture because the old source
does not automatically complete arbitrary selected supports. This does not
backport production code or change the selected variational space, but it
does not test old default admission behavior either.

Each native timing covers action construction, one H/S action, accepted
gradient preparation, one orbital direction, and its cached reverse. It is
**not a molecular SCF convergence time**. Every entry runs once per variant;
the current-only instrumented run is for attribution, not timing claims.
Numerical comparison includes all entries of 10 output/adjoint arrays, not
only norms or checksums. The current source's unrelated RI changes are not
used by this packed-exact fixture.

Release flags: GCC 14.1.0, `-O3 -DNDEBUG -g1 -fno-omit-frame-pointer`,
Eigen BLAS/LAPACKE enabled; 8 OpenMP cores in one LLC domain, BLAS single
threaded, identical affinity within the same Slurm allocation.

## Matched GitHub result

Job 1262586 completed successfully on anode50. Both variants used CPUs
80--87; all entries of the ten numerical blocks agree, with maximum
block-relative error 3.74e-15 and maximum absolute entry error 2.67e-15.
GNU time process RSS includes the fixture and all library workspaces.

| Structures | GitHub total (s) | Current total (s) | Speedup | GitHub RSS (MiB) | Current RSS (MiB) |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.6771 | 0.3943 | 1.72 | 160.7 | 58.7 |
| 3 | 0.7685 | 0.4112 | 1.87 | 160.3 | 75.5 |
| 20 | 0.7380 | 0.5020 | 1.47 | 160.9 | 75.6 |
| 132 | 0.7537 | 0.5057 | 1.49 | 161.5 | 76.4 |

| Structures | H/S action old/new (s) | Accepted gradient old/new (s) | Orbital direction + reverse old/new (s) |
|---:|---:|---:|---:|
| 1 | 0.0903 / 0.0658 | 0.2301 / 0.0720 | 0.3538 / 0.2481 |
| 3 | 0.0996 / 0.0650 | 0.2730 / 0.1430 | 0.3886 / 0.1882 |
| 20 | 0.0926 / 0.1370 | 0.2608 / 0.1499 | 0.3816 / 0.2117 |
| 132 | 0.0922 / 0.1242 | 0.2614 / 0.1499 | 0.3800 / 0.2278 |

This is a net improvement, not uniform improvement of every stage. The
20/132-structure H/S action still costs more than the old resident dense
carrier implementation. Generated-path planning and output projection are
the next targets there; the accepted and directional derivative chains
already more than offset that cost in this fixture. RSS is process peak,
not a claim that all metadata and consumers fit within the sigma slab limit.

## Isolated effect of this round

The GitHub comparison includes preceding low-memory changes, not just this
fusion. To isolate the new fusion, reuse saved predecessor job 1262522
(`48df68ab`) and compare with current job 1262574 on the same anode26,
8 LLC-local cores. No predecessor executable was rerun.

| Structures | Previous direction + reverse (s) | Fused direction + reverse (s) | Time reduction |
|---:|---:|---:|---:|
| 1 | 0.2866 | 0.2893 | -1.0% |
| 3 | 0.2355 | 0.1902 | 19.2% |
| 20 | 0.2794 | 0.2296 | 17.8% |
| 132 | 0.2795 | 0.2321 | 17.0% |

The one-structure case shows no measured gain. All measurements are single
runs; no significance claim is made for a one-percent difference. The
latest change to the panel-extent bound protects very unequal spin spaces
and leaves the CAS12 extent unchanged; final-source CTests are job 1262578.

The machine-readable matched results are `joint_direction_github12.csv`.
`joint_direction_evidence.tar.gz` preserves raw matched outputs, all ten
numerical arrays, build flags, binary hashes, affinities, regression and
independent-oracle logs, and the production source patch for this round.
The frozen executables remain in the Slurm artifacts on Hanhai22.
