# Exchange-reduced orbital generator contraction

## Scope

This identity applies to identical alpha and beta spin-string spaces. It does
not assume a complete VB structure space, an invariant selected space under
orbital transformations, or a negligible CI residual. The existing exact
structure-incidence test must establish one exchange parity for every selected
structure. Unequal spin spaces and mixed parities retain both spin contractions.

Let $P$ exchange the two string indices, so that $P X=X^T$. For two complete
carriers $L,R$ with the same parity $s=\pm1$,

$$
P L=sL,\qquad P R=sR,\qquad
P E^\alpha_{pq}P=E^\beta_{pq}.
$$

The orbital-generator cotangent is

$$
D_{pq}(L,R)=\langle L,(E^\alpha_{pq}+E^\beta_{pq})R\rangle.
$$

Orthogonality and involutivity of $P$ give the exact identity

$$
\langle L,E^\beta_{pq}R\rangle
=\langle PL,E^\alpha_{pq}PR\rangle
=s^2\langle L,E^\alpha_{pq}R\rangle,
\qquad D_{pq}(L,R)=2\langle L,E^\alpha_{pq}R\rangle.
$$

It also holds for a linear combination of cotangents, requested orbital-index
transposes, and the final upper-triangular Cholesky pullback. Spin-free sigma,
spin-free integral variations, and a common alpha/beta exterior transform and
its directional derivative preserve the parity. Therefore accepted, orbital
directional, coefficient-response, and residual carriers satisfy the premise.

## Bounded-panel realization

For disjoint beta panels $J$, compute

$$
D_{pq}(L,R)=2\sum_J\langle L_{:,J},(E^\alpha_{pq}R)_{:,J}\rangle.
$$

An individual panel generally has no exchange parity. Its reduced contribution
need not equal the corresponding two-spin panel contribution. Equality holds
after summing all panels; linear partial-image eviction is allowed because the
sum of the emitted left images reconstructs the full carrier.

The implementation omits beta spectator metadata, column gathers, and link
contractions in `generator_adjoint_streamed`. The tiled sigma uses weighted
exchange representatives internally; its emitted image is not automatically a
parity carrier. For complete resident output panels, first apply
$P_s=(I+sP)/2$ in place to all image components. Since the spin-free generator
commutes with $P_s$ and the right carrier has parity $s$,

$$
\langle L,(E^\alpha+E^\beta)R\rangle
=\langle P_sL,(E^\alpha+E^\beta)R\rangle.
$$

Each unordered product pair is updated once. Compact support is invariant under
exchange, so transpose entries are already present and no numerical array is
added. Evicted partial outputs retain the general two-spin contraction; they
must not individually use the parity identity. Retained accepted residuals use
the projected image. Alpha panels and small orbital
reductions retain their existing bounded storage. There is no new full carrier
or per-thread coefficient array. The alpha and beta traversal costs were
additive; this identity removes the beta traversal, not the asymptotic product
dimension or the physical structure response.

Accepted compact support indices are also immutable at a fixed orbital point.
Directional input, output, and projection tiles borrow the accepted CSR pattern
when their panel extent agrees. Missing retained data still use the same exact
support generator under the unchanged workspace budget.

## Verification

The independent oracle evaluates both spins with
`one_body_generator_adjoint`. Tests sum incomplete reduced panels for both
parities, multiple right arguments, requested transposes, and symmetric
pullbacks. A throwing beta-column reader proves the reduced contraction never
requests beta spectators. Existing unequal-spin and transient-reader tests
exercise the general two-spin contraction. Physical directional and response
oracles and a converged CAS(12,12) calculation qualify the integrated path.
