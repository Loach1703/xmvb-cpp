# Tiny spinor reference fixtures

`test_spinor_strings.cpp` generates these fixtures deterministically; no input
file or pseudorandom library is involved (random seed: none). The sine/cosine
formulas in the test source fully specify every coefficient. The reference implementation lives in
`tests/vbscf/unit/spinor_reference.hpp` and is test-only.

The complete `(2o,2e)` model has four spinors, two A and two B, and six
determinants. The `(6 spinors, 3 electrons)` model has three A and three B
spinors and 20 determinants; its two-body matrix contains transfers between
sectors whose A occupations differ by two. Vacuum and full-occupation cases
exercise the endpoints. A separate four-spinor, one-electron Hamiltonian obeys
time reversal and checks Kramers degeneracy. Generic random-looking complex
Hermitian matrices in the other fixtures are **not** claimed to obey time
reversal.

Spinor indices start at zero. Each reference determinant is an ascending list
of occupied indices, ordered lexicographically. Production masks are used only
to permute the completed dense reference for comparison. Creation and
annihilation signs come from insertion and deletion positions in the list.
No production string link or parity helper is called.

The one-body term is `sum h[p,q] a_p^+ a_q`. The exact two-body term is
`(1/2) sum g[p,q,r,s] a_p^+ a_r^+ a_s a_q`, with flattened tensor index
`(((p*n+q)*n+r)*n+s)`. The rightmost operator acts first. RI factors are
deterministic Hermitian complex matrices, and the *same* exact tensor used by
the dense oracle is reconstructed as `g[p,q,r,s] = sum_P L_P[p,q] L_P[r,s]`.
This gives pair exchange and Hermitian reversal, without imposing the invalid
independent reversal of either complex density pair. The diagonal identity
factor has the analytic two-body limit `N(N-1)/2`, including zero for `N=0,1`.

The helper `wedge_matrix(R, masks)` evaluates complex minors directly and
includes the vacuum minor equal to one. It is available for the selected
structure projection tests. A nonorthogonal fixture can use an upper
triangular complex `R` with positive real diagonal, `G = R^+ R`, and a sparse
three-column structure map `B`. Then `T = wedge_matrix(R, masks)` and the
independent dense carrier give `M = (TB)^+ (TB)` and `H = (TB)^+ H_orth (TB)`.
An occupied-pair overlap can be zero even when `R` is invertible; this tests
that the frame is judged by its full rank, not by one subset overlap.

All fixtures are intentionally tiny. They validate conventions and algebra,
not molecular X2C input, optimization, or large-carrier performance. Numerical
acceptance on Hanhai25 is recorded separately from this source description.

## Implemented projection checks

`test_spinor_projection.cpp` checks complex Gram whitening, reconstruction of
the original coefficients, shear-factorized exterior action and its adjoint,
complete and three-column selected spaces, exact/RI agreement, and equality
of projected bilinear forms with their full exterior images. It also checks
the independently evaluated compound Gram overlap, zero occupied-subset
overlap, singular-frame rejection, and workspace-budget rejection. The
reference tolerances are relative norm bounds of `1e-11` (strings) and `2e-11`
(projection); these are test tolerances, not molecular convergence criteria.

The production action stores no carrier-pair matrix or dense structure H/S.
It still uses complete fixed-N carrier amplitudes. The row cap is checked
before string allocation; the workspace budget covers amplitude blocks and
matrix-product scratch, not static input storage or allocator overhead.
Exact integrals remain an ordered `n_spinor^4` input tensor. This implementation
is not a solution to exponential carrier growth.

## Cluster acceptance — 2026-10-09

Hanhai25 Slurm jobs `278428` (clean build) and `278430` (final projection
test rebuild) completed successfully with the repository Release CMake
configuration and 16 compilation workers. All numerical tests ran on compute
nodes, with one OpenMP/BLAS thread: `spinor_strings`, `spinor_projection`, and
the existing scalar `orthogonal_direct_ci` regression all passed. Final CTest
elapsed time was 0.06 s; this is tiny-test runtime, not a molecular benchmark.
The fixture source is deterministic and contains no external downloaded data.
