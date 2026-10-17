# Exchange weights without mutable channel slabs

## Exact contraction

Let $A,B\in\mathbb R^{N\times P}$ be the response and accepted-state
one-body channel matrices for a generated tile. Each row denotes one
spin-string-product representative and one RHS. Here $P$ is the packed
orbital-pair channel count, not a number of VB structures. A fixed exchange
parity allows one representative of each alpha/beta exchange orbit. Its
multiplicity is $\mu_r=2$ for $a_r\ne b_r$ and $\mu_r=1$ for $a_r=b_r$.

The required integral-adjoint Gram matrix is

$$
G=A^T D B,\qquad D=\operatorname{diag}(\mu_r).
$$

Let $J$ select precisely the diagonal rows $a_r=b_r$. Since
$D=2I-J^TJ$, one has

$$
G=2A^TB-(JA)^T(JB)
 =2A^TB-\sum_{r:\,a_r=b_r}A_r^TB_r.
$$

Only the diagonal-row indices are needed. They are obtained during existing
target compaction, including generated targets outside the original support.
A target's own coefficient may be zero while its channel image is nonzero;
the original-support label must not be used to omit the correction. Arbitrary
selected VB structure spaces require neither closure nor a complete CAS for
this identity. Without exchange reduction use $G=A^TB$ unchanged. Odd exchange
parity excludes diagonal representatives, so its correction is empty.

## Row-partitioned and output-partitioned forms

For disjoint row partitions $R_t$,

$$
G=\sum_t\left(2A_{R_t}^TB_{R_t}
 -\sum_{r\in R_t:\,a_r=b_r}A_r^TB_r\right).
$$

The correction is applied to each existing small partial Gram before its
unchanged reduction. When the workspace instead partitions output channels,
the same identity applies to each disjoint block of $G$. No new channel slab,
per-worker numerical owner, full product tensor or selected-space matrix is
introduced. The shared diagonal-index list has at most one entry per target
beta group (before repeating it across RHSs), hence O(groups) metadata.

## Directional response remains unchanged

The unweighted accepted channel slab is also needed for a directional
Hamiltonian image. Previously it was weighted for $G$, then unweighted before
the direction contraction, or its directional image was unweighted afterward.
The new representation keeps $B$ unweighted throughout all temporal requests.

$$
W=A K^T+B(\delta K)^T,\qquad W_{\mathrm{force}}=B(\delta K)^T.
$$

The existing exchange reconstruction of these output images is unchanged.
Thus no CI response, orbital forcing, integral adjoint or normalization term is
discarded. The change only relocates multiplicity handling from large channel
slabs to a small Gram plus exact diagonal-row corrections.

## Cost and numerical qualification

The dominant Gram multiplication remains O(N P^2). Removing slab weighting
and unweighting eliminates O(N P) passes; the diagonal correction introduces
O(d P^2) work for d diagonal rows. This is not an asymptotic reduction of the
complete sigma/HVP algorithm and need not be faster for small, diagonal-heavy
tiles. No empirical threshold selects rows or changes accuracy.

The two expressions have different floating-point summation order. Retention
requires the existing independent dense/selected action, adjoint and HVP
finite-difference tests at their original tolerances, followed by a converged
C12 end-to-end test. Mathematical equivalence alone is not a timing result.
