# Generated direct-CI sigma/VJP and the triangular-support certificate

Parent checkpoint: `fb4b0bc4`. Numerical runs use Hanhai22/anode26, 8 bound
OpenMP cores, BLAS 1, unchanged Release flags and physical tolerances.
Historical executables are not rerun. This is a staged streaming migration,
not a completed removal of broad-support derivative carriers.

## Shared generated channel graph

`DirectCiSigmaAction::apply_tiled_with_adjoint()` uses the existing generated
paths and gather graph. It forms the integral cotangent from the same channel
scatter used by sigma. Equal inputs share their scatter; distinct inputs
reuse the consumed input halo and future K-times-channel slab. The Gram is
consumed before that slab is overwritten. One shared pair-kernel cotangent
has disjoint worker-owned rows; there is no per-worker spin/product array.
The numerical slab target remains 64 MiB. The returned active-integral
cotangent and path metadata are separate, as in the established dense order.

The API accepts generated readers, not only resident product arrays. The
sigma-only public API delegates to this same implementation. It preserves
the full output outside the input support unless the caller explicitly
requests a certified exchange projection. Both equal-parity sectors and
unequal spin spaces are tested; the integral cotangent weights exactly one
Gram operand by each exchange-orbit multiplicity.

## Bounded generated-support probe1261518

Constructed exact 12-AO Hamiltonian, 924 strings/spin. Input amplitudes are
generated analytically on the complete-covalent downset. This is an algebraic
kernel stress test, not a physical converged VBSCF calculation. The native
process contains neither a complete input product array nor complete sigma;
its output consumer contracts partial sigma columns to a scalar.

|Quantity|Native generated run|
|---|---:|
|Generated channel targets|608751|
|Complete product targets|853776|
|Batches|32|
|Peak counted slab payload / bytes|34899120|
|Preparation / s|0.011227|
|Analytical reader / s|0.083102|
|Channel computation / s|0.144446|
|Scalar output consumer / s|0.274992|
|Complete generated sweep / s|0.517557|
|Process peak RSS / MiB|42.7813|

A separate process checks the same new code against the complete-carrier
algebraic oracle: energy error6.61e-14, one-electron cotangent error2.69e-16,
pair-kernel cotangent error1.55e-15. Its dense sigma/adjoint subtotal is
0.233600 s; process RSS is87.1719 MiB, including its generated sweep and
explicit reference arrays. These are different reader/output contracts,
not a matched total-time or memory speedup comparison. In particular the
simple serial scalar consumer dominates the native probe. This primitive
is correctness/resource-qualified, **not** a universal speedup result.

Values and stage counters are parsed into `tile_adjoint_kernel.csv`; no
sampled CPU percentages or allocation-profiler claims are inferred.

## Discarded resident-carrier adapter

An experimental adapter fed existing resident accepted coefficient arrays
through the generated interface. It did not remove those arrays. Job1261520
passes11 tests and9 complete240 runs but gives no meaningful end-to-end gain.
FeCl2 job1261524 takes16 steps/269 HVP,12.251995 s,404.996 MiB, final energy
-2181.624464847211 Eh. The parent control1261472 took13 steps/235 HVP and
10.241470 s. This adapter is deleted, not kept as a switch or fallback.
Current production preparation continues to use its established exact
factor/dense contraction order until its readers and residual consumers
are migrated together. The rejection concerns that premature integration,
not the algebraic correctness of the new generated primitive.

## Triangular support, rather than assumed structure-space closure

The Cholesky factor and its relative direction `L=delta-R R^-1` are upper
triangular. Upper occupation/replacement generators preserve the transformed
downset S. Both selected-coefficient and orbital-representation tangents
therefore remain in S; Hamiltonian channel intermediates use its one-generator
closure I. Differentiating this triangular transform does not itself require
a larger second closure.

The Cholesky overlap pullback uses only the upper entries of the oriented
residual-generator cotangent. These entries and their directional derivatives
are unchanged if the residual is projected to S. The full oriented cotangent
is generally changed: this is not a license to truncate its generic interface,
nor to project residuals into the selected VB structure basis.

Tests cover4/5/6 active orbitals, covalent and arbitrary ionic/mixed support,
multiple RHS and unequal spins. The upper-generator tangent has exactly zero
support leakage. The discarded full-generator norms are0.30048,1.41792,
1.58417 and0.25751, confirming that nontrivial components really are omitted,
while physical overlap pullbacks agree within3e-11. Job1261552 also checks
the differentiated physical overlap pullback against the complete-residual
oracle. There is no change to production residual storage yet.

The two proofs and their cost/storage scope are documented in subsections4.2
and4.3 of `article/source_column_hamiltonian.tex`; the rebuilt PDF was inspected.

## Final code verification, after deleting the adapter

Job1261533 passes11/11 focused tests and9 complete240 optimizations.
Job1261552 passes11/11 tests, including the final added triangular-support
derivative and differentiated physical overlap-pullback comparisons.
No optimizer policy, physical convergence threshold or coordinate update is
modified by this stage.

|240 space|Steps|HVP|Final energy / Eh|SCF / s|RSS / MiB|
|---|---:|---:|---:|---:|---:|
|COV|3|20|-343.239211501260|6.465714|1514.43|
|ION(0:1)|3|27|-343.407736636065|8.655557|1519.70|
|ION(0:2)|4|32|-343.443415326852|10.830612|1521.13|
|ION(0:3)|4|32|-343.446256564846|11.121905|1527.10|
|FULL|4|30|-343.446302545201|9.425874|1520.87|

Steps/HVP counts agree with saved parent data, energies within3e-12 Eh.
Selected1/3-structure dense/Davidson cases all take3 steps/21 HVP.
Final unequal-spin FeCl2 control1261545 passes8 geometry/coordinate tests:
15 steps/249 HVP, -2181.624464861484 Eh,11.185807 s SCF,428.922 MiB RSS.
Energy agrees with the parent within3.0e-9 Eh, but its action count is higher
than the saved13-step/235-HVP parent and12-step/211-HVP correct-metric result.
The source of that trajectory difference is not established here; it is
not attributed to rounding or excused by a precision change.

All are exploratory single-run times on the same node with other allocations
present; iteration/action counts and process native RSS are recorded separately.
The complete source-state CAS16/16 HVP speedup of the preceding commit
(100.29 to76.26 s) is retained as saved evidence, not rerun or relabeled as
a new complete-SCF speedup.

## Next production migration

Use generated coefficient readers in the accepted product consumer and
consume the sigma/residual on S before exterior/structure backprojection.
Accumulate the certified upper generator VJP there alongside the integral
cotangent; leave generic full-oriented generator APIs intact. Reuse this
same graph for orbital/coefficient directions, then delete retained broad
accepted/directional product arrays. Do not merely wrap resident arrays or
discard residuals by assuming selected-space closure.

Raw evidence, including the rejected integration's traces, is archived in
`tile_adjoint_evidence.tar.gz`; final/candidate CSV files are named separately.
Development failures1261510/1261514/1261517 are retained too: a test identifier
collision and a task-copy error in the subsequently deleted dense adapter were
caught before qualification. They are not final-code test passes. Duplicate
executables and unrelated metric-point binary matrices are excluded from the
archive; their source/binary hashes and text diagnostics remain.
