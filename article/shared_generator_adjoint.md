# Shared-right generator adjoints

For a left beta panel $L_j$ and generated right arguments $R_r$, the oriented
orbital moment is

$$
D_{pq}(L_j,R_r)=\sum_{a,b,a',b'}
 (L_j)_{ab}\,\langle ab|E_{pq}|a'b'\rangle\,(R_r)_{a'b'}.
$$

The sum over $b$ is restricted to the supplied panel. Alpha links keep $b=b'$;
beta links keep $a=a'$. This is a spin-string link contraction, not an explicit
determinant-pair matrix. A requested orientation and scale act linearly on the
result. Define the same right functional for several left panels,

$$
G_j=\sum_r s_r\,\mathcal O_r D(L_j,R_r),
\qquad \mathcal O_r D\in\{D,D^T\}.
$$

The $G_j$ must remain separate. Their right beta-column reads, signed-link
topology and spectator gather addresses are identical and need be generated
only once. Dot products with different $L_j$ remain necessary: linearity does
not permit returning one sum as both responses. Each worker owns only one
orbital-square partial per left panel; all carrier inputs remain borrowed.

For the physical coupled product jet, the two left panels are the orbital
Hamiltonian direction and the CI Hamiltonian image. Both contract against the
same accepted $Y$. In an uncached accepted point, the orbital panel is
$H\dot Y+\dot H Y$, while a separate $HY$ panel still has to contract with
$\dot Y$ and $Q$. These latter independent right arguments are not discarded.

For a resident accepted residual $R=(H-E)Y$, the tail moments are

$$
G_p=-E D(\dot Y,Y)+D(\dot Y,R)^T,\qquad
G_q=-E D(Q,Y)+D(Q,R)^T.
$$

Thus both use the same two right arguments and can share another traversal.
The existing selected-state gauge correction, upper-triangular orbital
pullback and subsequent AO adjoint remain unchanged.

Exchange reduction is applied only under its existing exact equal-parity
certificate. It replaces the full alpha-plus-beta sum by twice the alpha sum
after all left panels have been accumulated. Neither an isolated panel nor
unequal spin spaces are assumed to have this property. Symmetrization and
upper-triangle extraction are performed independently on each $G_j$.

If the right panels occupy $B_R$ bytes, a batch of $m$ left panels uses
$B_R+O(mtn_{\mathrm{act}}^2)$ temporary storage, where $t$ is the worker count.
The left panels are non-owning views. Right reads and link discovery lose a
factor $m$, but arithmetic with the distinct left arguments does not. Therefore
this is a reuse optimization, not a claim that the full contraction's formal
polynomial order has decreased.

The single-left interface delegates to the same kernel; no duplicate or
compatibility contraction engine is retained. Product-jet output panels also
borrow their existing Hamiltonian blocks instead of copying them before these
read-only contractions. Hanhai22 job1264759 passed the independent multi-left
oracles and full coupling regression, with unchanged four-step/36-HVP convergence.
The scoped generator calls decreased704 to416, but native time changed only
29.0093 to28.8680 seconds and instrumented total time increased slightly. The
reuse is real; it is not evidence of a substantial end-to-end speedup. The
larger finite-transport pullback remains the next measured target. Raw evidence
is under `results/diagnostics/c12-shared-generator-20261017/`.
