# Fused finite-transport generator moments

## Required derivatives

The finite CI transport gradient contains a target overlap, a connection term,
both derivatives of the current metric, and an origin-chart term. None may be
discarded by a small-residual or step-size test. The factor representation
expresses a target image as a sum of outer products of single-spin panels.
For one spin, define an upper one-body moment and its chart pullback by

$$
J_{pq}(A,B)=\sum_{i,j} A_{ij}(E_{pq}B)_{ij},\qquad
P_R(J)=\operatorname{upper}(J R^{-T}).
$$

Here $A$ is the seed panel, $B$ the image panel, $E_{pq}$ the single-spin
one-body generator, and $R$ the upper orbital transform. Occupation and signed
single-excitation links evaluate $J$ without forming a compound matrix or a
spin-product carrier. Only the entries that can contribute to the upper chart
are generated, as in the established single-term implementation.

The pullback is linear in its moment, so for fixed $R$,

$$
\sum_\ell s_\ell P_R[J(A_\ell,B_\ell)]
=P_R\!\left[\sum_\ell s_\ell J(A_\ell,B_\ell)\right].
$$

Thus each shear link can produce the sum of its required bilinear dot products
in one visit. Their arithmetic remains necessary, but rediscovering the same
occupation/shear topology and solving the same upper-chart pullback repeatedly
is unnecessary.

## Factor-panel application

For the alpha current chart, the four terms already present in the finite
transport implementation are

$$
J(U_b,Z_a^{(1)})+J(U,Z_a^{(2)})
-J(U_c,Z_a^{(m)})-J(Z_a^{(m)},U_c).
$$

The beta current chart analogously requires

$$
J(V,Z_b^{(1)})+J(\dot V,Z_b^{(2)})
-J(V_c,Z_b^{(m)})-J(Z_b^{(m)},V_c).
$$

$U_b,V$ are the first target factors, $U,\dot V$ the second target factors,
and $U_c,V_c$ the transported current factors. $Z^{(1)},Z^{(2)},Z^{(m)}$ are
their existing seed contractions against the selected left coefficients.
The two negative metric terms have opposite operand orders; they are not
interchangeable and neither is removed. The origin-chart terms remain
$P_{R_o}[J(Z_a^{(1)},U)]$ and $P_{R_o}[J(Z_b^{(2)},V)]$.
The original notation corresponds respectively to `ub,u,uc,v,dv,vc` and
`first,second,metric` in `transport_tiles.cpp`.

Per factor panel, ten generator pullbacks become four. Different spin
topologies remain distinct; equal alpha/beta dimensions are not assumed.
The transpose image still consumes $Z_a^{(1)}$ only after every chart consumer
has finished reading it.

## Directional and parallel consistency

The common kernel also retains the directional rule

$$
\dot J(A,B)=J(\dot A,B)+J(A,\dot B),\qquad
\dot P_R=\operatorname{upper}\!\left[
\dot J R^{-T}-J R^{-T}\dot R^T R^{-T}\right].
$$

Each off-diagonal upper generator corresponds to one immutable shear list.
Assigning one owner to each such entry needs no atomics and no thread-private
carrier or orbital-square reduction. The occupation diagonal computes each
row dot once, accumulating only an orbital-length vector per worker. Both
directional and nondirectional calls share this implementation, including
cooperative execution inside an existing OpenMP team.

The existing factor budget is $16w(U_\alpha+U_\beta)$ doubles. With all three
seed pairs retained, factor panels peak at at most $7w(U_\alpha+U_\beta)$
doubles during seed construction, leaving the established conservative budget
unchanged. Additional shared moments are orbital-square, and diagonal worker
reductions are $O(tn_{\mathrm{act}})$ doubles. No string-product or string-square
allocation is introduced.

These identities reduce repeated work and permit exclusive-output parallelism;
they do not change the formal order of the necessary bilinear arithmetic.
Qualification must compare the selected-space transport against its independent
cross-overlap derivative oracle, directional pullbacks against finite
differences, then full C12 convergence against saved old measurements.

## Measured qualification

Hanhai22 job1264773 passed all nine regression targets, including selected
cross-overlap/minor oracles, unequal-spin transport and directional finite
differences of fused moments. Both ordinary and existing-team OpenMP execution
are covered. The native CAS(12,12) LOFLEA task converged in four steps with
36 HVPs to $-399.106646267328$ hartree. Application wall time was27.623396s
and GNU-time peak RSS2679.746MiB, versus the saved immediate predecessor's
28.868012s and2682.656MiB. No old executable was rerun.

In the separate phase run, generator pullbacks fell440 to176 calls and
2.317848 to0.190112s; the full transport-gradient stage fell2.685966 to0.593678s.
The reduction combines moment fusion and exclusive-output parallel execution;
it does not isolate the contribution of either change. Total instrumented
application time fell29.068800 to26.449473s. Native and instrumented timings
are separate measurements, and nested scopes must not be added together.
The per-step HVP counts remain6,17,10,3. Thresholds, state space, coordinates,
physical coupling and model policy are unchanged. Raw evidence is stored under
`results/diagnostics/c12-transport-moments-20261017/`.
