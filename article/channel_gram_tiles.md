# Two-dimensional channel Gram contraction

In one generated spin-string tile, let $A,B\in\mathbb R^{m\times d}$ be the
response and accepted one-body channels. Spin-exchange multiplicities, if
applicable, have already been applied to $B$. The exact kernel cotangent is

$$
G\mathrel{+}=A^T B.
$$

Every output entry is required. No CI response or density component can be
discarded. However, row-stripe parallelism makes each of $P$ workers multiply
only roughly $d/P$ output rows against all $d$ columns of $B$. Its logical input
volume is approximately $(P+1)md$ doubles: the left columns are partitioned,
but the whole right panel is requested by every worker.

Partition both output axes into $b$ blocks. A block $(I,J)$ is exactly

$$
G_{I,J}\mathrel{+}=A_{:,I}^T B_{:,J}.
$$

Each block has one owner and only reads its own two channel subsets. The
logical input volume becomes approximately $2bmd$ doubles, without changing
the $2md^2$ floating-point operation count. For $P=32$, an $8\times8$ grid
gives about half the row-stripe input volume, before cache reuse and BLAS
packing are considered. This is not a measured DRAM-traffic claim.

The block count per axis is the smallest power of two whose square covers
the actual OpenMP team. The extent follows from $d$; it is not a new tuned
orbital or convergence parameter. For power-of-two team sizes, the grid has
one or two output tasks per worker, except for dimension truncation.

Only Eigen views of the existing shared row-major $A,B$ slabs and the existing
small output cotangent are used. There is no new numerical workspace,
per-thread carrier, replicated integral array, or atomic reduction. The
existing end-of-Gram barrier still protects the subsequent reuse of channel
storage. General unequal-spin sectors and both exchange parities use the
same Gram identity. Accuracy, support generation, Hamiltonian application,
sigma reconstruction and all optimizer decisions are unchanged.

Opt-in timing records the whole-team Gram wall scope (nested in the existing
tile compute scope) and its output-block count. Native runs do not execute
the timing-only synchronization. Old versions are not rerun.
