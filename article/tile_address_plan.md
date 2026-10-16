# Relocating spin-string addresses once per tile

For a tile with $H$ alpha halo rows and $b$ right-hand sides, let $h(a)$ be
the halo position of alpha string $a$ and $c(\beta)$ the local source-column
position of beta string $\beta$. Column-major input packing gives

$$
j(a,\beta,v)=h(a)+H\bigl(bc(\beta)+v\bigr),\qquad 0\le v<b.
$$

For a target beta group $g$, every alpha target and RHS uses the same beta
replacement edges $\nu=(\beta_g\leftarrow\gamma_\nu)$. Relocate their
addresses once:

$$
o_{g\nu}=Hb\,c(\gamma_\nu),\qquad
j(a,\gamma_\nu,v)=\bigl(h(a)+Hv\bigr)+o_{g\nu}.
$$

The hot scatter then reads a shared compiled offset instead of loading a
source-string index, looking up its local column, and multiplying by both
RHS width and halo height for every alpha target. The owned beta column is
similarly relocated once per group. The response and accepted scatters reuse
this same address table; they still use the original single halo and separate
scatter passes, preserving the fast established GEMM/channel layout.

Address relocation work is $O(gd_\beta)$ per batch, rather than repeating it
for $O(mbd_\beta)$ channel visits, where $g$ is the number of beta groups,
$d_\beta$ their single-replacement degree and $m$ the generated target count.
The physical channel arithmetic still has its original scaling. Storage is
only the bounded $g\times d_\beta$ integer table, never a table on the full
alpha/beta product space. Each input panel is below 64 MiB, so offsets measured
in doubles fit safely in a signed 32-bit integer under the fixed budget.

A sorted distinct alpha halo with cardinality $U_\alpha$ is the identity
mapping $h(a)=a$. Scatter and transpose gather then use $a$ directly. This
does not assume a full selected structure space: it is a property of one
generated halo, independently for each batch. Partial halos retain their
exact maps. Diagonal support masks and all neighbor contributions are unchanged.

No electronic integrals, coefficient arrays, or thread-local numerical buffers
are added. Tile geometry, arithmetic precision, the physical HVP and structure
responses, and all optimization decisions are unchanged. A shared-scatter
two-halo candidate and a column-major image candidate were also investigated
but not retained, because their complete native runs did not improve on the
saved baseline. Faster-looking local operations alone do not justify a route.
