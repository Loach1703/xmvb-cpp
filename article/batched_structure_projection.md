# Independent-image batching in selected-structure projection

For a beta output panel $J$, let $Z_1$ and $Z_2$ be independent images, such
as the orbital-direction and coefficient-direction Hamiltonian images. Both
need the same linear selected-space projection. With $\mathcal E$ the
original sparse structure expansion, this projection is

$$
\mathcal P_R(Z)=\mathcal E^T\operatorname{vec}
\left(T_\alpha^T Z T_\beta\right).
$$

Apply the alpha adjoint to the concatenation $[Z_1\ Z_2]$ once. Its columns
are independent; the supported implementation repeats exactly the same beta
column pattern for each image. For original touched product $g=(a_g,b_g)$,
evaluate both projections in a single support traversal:

$$
v_{g\ell}=\sum_{j\in J}
(T_\alpha^T Z_\ell)_{a_gj}(T_\beta)_{jb_g},\qquad \ell=1,2.
$$

The two values retain independent selected-structure output columns. This
removes a repeated visit to original product indices, projection-tile lookup,
beta-row resolution and cooperative-loop launch, not either physical action
or a cross term. The derivative projection still keeps both
$\dot T_\alpha^T Z$ and $T_\alpha^T Z\dot T_\beta$. No selected-subspace
closure is assumed, and zero-amplitude products are not screened.

## Ownership instead of a second copy

The old const-owning-matrix argument first converted a block expression to an
owning temporary, then copied it again into the mutable transform workspace.
An owning by-value argument needs only the first copy. A caller's full local
image can instead be moved after its last generator/diagnostic consumer.
Accepted cache storage remains separate and immutable: moving an unpacked
image must never move or alter its compact cache owner.

For joint overlap images, load the two directions directly into adjacent
blocks of one shared panel, rather than allocating and concatenating two
independent panels. All oriented generator contractions finish before its
destructive projection. Diagnostic overlap and norm consumers run before this
final physical projection, so capture does not require a copy of both
direction panels. It changes consumer ordering, not the mathematical operator.

## Memory boundary and acceptance

In the optional resident projection, allocate a left-direction image only
for a target that actually receives a differentiated projection. A CI-only
or integral-forcing target has no such component: its final beta operation
is an ordinary adjoint, not a directional adjoint of an explicitly allocated
zero carrier. This omission follows from the requested operator components,
not from an amplitude, residual or available-memory test.

Accepted left-direction terms share the same projection as well. Unpack the
immutable accepted $[Y,(H-E)Y]$ tile once, replace its disposable second block
by $HY$, and directionally project $[Y,HY]$ to the independent overlap and
Hamiltonian outputs. This eliminates a residual copy, a separate reconstruction
temporary, the later dense $Y$ reread and one more support traversal. Diagnostic
representation forcing consumes its $HY$ copy before the combined panel is
moved. The original accepted cache remains untouched and still serves all
oriented generator contractions and physical response consumers.

The original panel extent and budgets remain. A joint Hamiltonian callback
has at most three output components and two copied projection components;
the existing six-alpha-column charge bounds this workspace. A cached joint
callback can move its two components without any new carrier. Overlap batching
uses the same two direction components already required, with no per-thread
carrier, additional accepted cache or full product-space mandate.

The formal arithmetic order is unchanged: alpha transforms and scalar
projection products for both images remain necessary. The proposed gain is
from shared traversal and ownership, not a reduction in Hamiltonian-action
count or an asymptotic exponent. Qualify independent outputs, zero directions,
unequal spin spaces, both exchange parities, nonclosed support, resident and
evicted tiles, diagnostic capture and finite-difference HVPs before evaluating
native end-to-end convergence on Hanhai22. Existing old results are not rerun.
