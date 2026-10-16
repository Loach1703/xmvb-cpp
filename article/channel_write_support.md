# Exact row-range restriction of channel transpose scatter

**Status:** this is an exact but rejected implementation experiment, not the
retained production algorithm. C12 measurements did not establish an
end-to-end speedup. The experimental branches and their metadata were removed;
only independent sparse-output regression tests and archived evidence remain.
See `results/diagnostics/c12-beta-range-20261017/README.md`.

For one generated beta target group $G$, let $A_G=\{a_i:i\in G\}$ contain
its actual alpha-string indices. The channel image $W_{ip}$ contributes to a
beta transpose edge $e:G\to b$ as

$$
y_{ab}\mathrel{+}=s_e\sum_{i\in G}\delta_{a,a_i} W_{i,p_e}.
$$

The beta one-body operator changes the beta string, not the alpha spectator
index. Therefore this contribution is identically zero unless $a\in A_G$.
In particular it vanishes outside the closed interval

$$
[\ell_G,u_G]=[\min A_G,\max A_G].
$$

The original output-row set $O_b$ is already sorted. Restricting its scan to
$O_b\cap[\ell_G,u_G]$ is exact for every channel value, response RHS, exchange
parity and selected structure space. Keep the original target-index check
inside this interval: a range bound must not turn holes into real targets.
No monotonic ordering of the caller's unique strings or generated targets
is assumed. Only the existing sorted halo/output indices are searched.

The spectator restriction is **not** valid for alpha transpose scatter,
which changes the alpha string. Its distinct valid restriction follows below.
Neither selected-space closure nor CI stationarity is used, and no Hamiltonian
or structure-response term is screened by amplitude or residual size.

## Alpha reverse one-body closure

Let $N_\alpha(a)$ contain the channel-source indices referenced by the single
links of output alpha string $a$. Its transpose support for target group $G$ is

$$
H_G=A_G\cup\{a:N_\alpha(a)\cap A_G\ne\varnothing\}.
$$

Outside $H_G$, neither a diagonal channel nor any alpha single-link channel
can contribute. For each target string $t$, precompute the inclusive bounds
of $\{t\}\cup\{a:t\in N_\alpha(a)\}$ by explicitly reversing the immutable
link table. The minimum/maximum across $t\in A_G$ then bounds $H_G$.
This requires no assumption about ordering or reciprocal adjacency. It skips
an entire output row's link queries, not just a beta spectator lookup. Keep
the original exact map checks within the bound, including all holes.

The two bounds per string are shared with the accepted spin connections and
their integral directions; they are not reconstructed per HVP. Group bounds
are accumulated during existing target compaction. The additional metadata
is $O(U_\alpha+G)$, never $O(U_\alpha U_\beta)$ or per-worker carrier storage.
One clipping implementation serves both the wider alpha closure interval and
the narrower beta spectator interval; it must not confuse their meanings.

## Work and memory

Compute $\ell_G,u_G$ while the existing target-compaction loop already reads
each generated index. The beta-only experiment added two integers per group;
the alpha experiment added two more per group and two per immutable spin
string. Neither needed a new target list, sorting pass or product-space map.
Empty groups have an empty
interval and contribute no beta updates. Dense identity-halo writers clip
integer bounds directly; compact sorted output rows or nonidentity halos use
constant-time disjoint-interval tests first, then lower/upper bounds only when
their current range partially overlaps a group boundary.
Worker ownership and channel/reduction order are unchanged.

Previously the beta contribution visited $|O_b|$ candidates for every edge;
the new visit count is $|O_b\cap[\ell_G,u_G]|$. Full broad groups may give no
saving, so the worst-case asymptotic order is unchanged. Gain depends on the
actual generated support, especially exchange-representative or split target
groups. A bounded diagnostic records candidates and surviving range visits,
not numerical nonzero amplitudes. Its statistics are archived; the counters
are removed from the retained path. Native convergence timing must remain
separate from the diagnostic run.

## Affine addressing experiment

If the actual group is the consecutive ordered list
$a_i=\ell_G+i$ for $0\leq i<|G|$, its target row is exactly

$$
t_G(a)=\begin{cases}
\operatorname{begin}(G)+a-\ell_G,&\ell_G\leq a\leq u_G,\\
-1,&\text{otherwise}.
\end{cases}
$$

An explicit group-by-halo index map is then redundant. Unordered groups that
contain every integer in the interval can be put in that order by in-place
cycle placement in $O(|G|)$ work. The structural-support label must move with
its index. No such simplification applies merely because a group's bounding
interval is wide: holes require an exact map.

On the measured support, natural consecutive groups were rare. Canonicalizing
complete intervals enlarged coverage only modestly. Neither variant established
an end-to-end gain, so the affine specialization was removed too. This
experiment changes no worst-case asymptotic order and is not a new contraction
scaling result.
