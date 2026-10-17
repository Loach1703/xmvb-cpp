# Lifetime-exact workspace bounds for streamed CI contractions

Status: implemented and numerically qualified on Hanhai22. The measured
contraction work is lower, but the two native runs do not establish a stable
end-to-end speedup over the saved baseline.

## Simultaneously live data

For a generated batch, let $r$ count one-body target products, $h$ the union
of required alpha strings, $b$ the union of beta input columns, $g$ the beta
target groups, $w$ the RHS width and $d$ the Hamiltonian channel count.
Each channel slab has $wrd$ doubles. The input halo has $whb$ doubles and
the target-address map has $hg$ integers.

The response/state readers execute in this order:

1. Read response coefficients into the shared halo.
2. Scatter them into the response channel slab; synchronize all workers.
3. Refill the **same** halo with the accepted state.
4. Scatter into the accepted channel slab; synchronize before the Gram.
5. Form the integral cotangent and Hamiltonian channel images; gather outputs.

No subsequent consumer needs the original response halo. Thus the memory
upper bound must charge one input halo, not two. With $s=2$ channel slabs
for an ordinary request and $s=3$ for an integral direction, the bounded
numerical slab and target map require

$$
8w(hb+srd)+4hg\leq B,\qquad B=64\times 2^{20}\ \mathrm{bytes}.
$$

The maximum row count before accounting for the halo is therefore bounded by
$\lfloor B/(8swd)\rfloor$, not half that number. The existing batch planner
then tightens this bound using the generated beta union, the conservative
alpha-neighbor union and the target map. A split column resumes its original
support cursor; it is not regenerated or screened by coefficient size.
The existing runtime check of the actual directional payload is unchanged.

This changes neither the allocated 64 MiB slab capacity nor the individual
reader/output cache budgets. It removes a nonexistent concurrent allocation
from the planning proof. Fewer batch boundaries can avoid repeated halo
reads, address construction and synchronization without any new persistent
state or thread-private coefficient array.

## Why the two-slab direction prototype was rejected

An alternative overwrote directional channels using borrowed micro-GEMM
scratch in the dead halo. It passed independent numerical oracles and ten
regression tests, including unequal spins and tiny halos. On the matched
C12 convergence run it reduced batches but did not beat the saved end-to-end
baseline, even after fusing two GEMMs into a single copyback. It is not
retained. The lifetime-exact bound above keeps the established out-of-place
GEMM and row-contiguous output, rather than paying for in-place copyback.

The bound covers the named shared numerical slab and target map, not all
application allocations. Process peak RSS must still be measured separately.

## Matched C12 measurements

Use the original LOFLEA CAS(12,12), 132 covalent structures, HAO, exact
Libcint/cc-pVDZ, CIAH and Davidson. On Hanhai22 anode40, all native cases use
32 physical cores with affinity 60-91 and one BLAS thread. The saved job
1264636 is not rerun. New jobs 1264676 and 1264685 use the same new binary.

| Native metric | Saved | New | New binary repeat |
|---|---:|---:|---:|
| End-to-end / s | 30.863387 | 31.677382 | 30.361555 |
| SCF / s | 28.207024 | 28.794006 | 27.658785 |
| Peak process RSS / MiB | 2824.676 | 2730.438 | 2775.211 |
| Accepted steps | 4 | 4 | 4 |
| Physical HVP / CI responses | 36 / 36 | 36 / 36 | 36 / 36 |

Final energies agree within 2e-12 Eh; unchanged finite-difference HVP error
is 5.81e-11 and independent structure-response image/adjoint errors are
3.40e-16/1.63e-15. Nine regression tests pass. No precision or optimizer
policy changes are involved. The only post-test edit to the added unit cases
clarifies a comment after removal of the scratch prototype.

Separate phase-enabled runs measure 3109 -> 1639 batches and
1,964,120,314 -> 1,420,939,020 logical source elements. Halo-read time is
3.412612 -> 2.644525 s; coupled HVP/structure scope 12.414762 -> 11.990040 s.
These overlapping scopes are attribution, not independent speedups to add.
Actual tile payload is 35.395 -> 60.944 MiB, still below the unchanged
64 MiB slab capacity; process RSS is not this payload.

The simpler lifetime bound is retained. It proves fewer redundant reads but
does not support claiming a large wall-time gain. Independent orbital and CI
response sweeps still repeat the accepted channel scatter. Their temporal
fusion, preserving two cotangents and two coarse-model images, remains the
next unimplemented structural opportunity.
