# Exact two-sided channel projection: theory and rejected storage experiment

The identity below is exact. Its compressed-output prototype passed complete
sigma, directional response and integral-adjoint tests, but is not retained in
production: the reduced multiplication work did not yield a sufficiently
clear end-to-end advantage over the simpler executor to retain its additional
ragged image addressing. The production executor retains the earlier
bounded full-channel output tiles and the independently useful support/closure
DAG schedule. There is no extra output mode or compatibility implementation.
The allocation/work formulas below describe the tested proposal, not current
production storage.

Let $F=TW$ embed the selected structure coefficients in the orthogonal
spin-product carrier. Let $P_S$ select its exact symbolic support, so that
$P_SF=F$. This support projector is not a projector onto the selected VB
linear subspace. Write the orthogonal Hamiltonian action in its existing
one-generator form,

$$
H=E^T K E,\qquad Y=P_S Y.
$$

Here $E$ denotes the spin-summed one-body scatter, including the existing
one-electron correction in $K$. At a generated intermediate product $t$,
the existing certified channel selector $J_t$ satisfies

$$
(E P_S)_t=J_t J_t^T(E P_S)_t.
$$

The complete-image algorithm computes $KJ_t a_t$, where
$a_t=J_t^T(EY)_t$. If the final consumer needs only the selected image, its
exact contraction is instead

$$
P_S H Y
=\sum_t (E P_S)_t^T J_t
\underbrace{(J_t^T K J_t)}_{K_t}a_t.
$$

Thus the full $P$-channel image is unnecessary: scatter, multiplication and
gather can all retain the same $k_t$ channels. A target with no restriction
has $J_t=I$ and uses the ordinary kernel. A certified frontier envelope uses
its exact selector on both sides. The proof concerns the output projection;
it does not assert closure of selected VB structures under $H$ or $T$.
An API consumer requesting the complete $HY$ still requires the rectangular
$KJ_t$ product and every output channel. Values outside $S$ in the selected
output mode are unspecified and must not be used as a complete sigma image.
The actual selected-structure action is $F^T(P_SH P_S)F=F^THF$;
Davidson coefficients and residuals remain in the original structure basis.

The same identity applies separately to $H\delta Y$ and $\delta H Y$.
The integral cotangent already has both selectors,

$$
\overline K=\sum_t J_t a_t b_t^T J_t^T,
$$

so neither its value nor the physical structure response changes. No
gradient or Hessian tolerance is involved. For $m_g$ rows sharing a selector
of width $k_g$, image multiplication changes from
$O(\sum_g m_g P k_g)$ to $O(\sum_g m_g k_g^2)$, and the output channel slab
changes from $\sum_g m_g P$ to $\sum_g m_g k_g$ doubles. This is an exact
support-dependent work reduction, not a universal lower exponent: groups
with $k_g=P$ retain their original cost. The existing bounded tile planner
charges all live ragged slabs and selector kernels; no product-sized cache
or per-worker numerical slab is introduced.
