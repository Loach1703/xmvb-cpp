# Arbitrary selected-space CI inverse-model audit

2026-10-17. C12/LOFLEA CAS(12,12), HAO, exact Libcint/cc-pVDZ,
ISCF=9/CIAH and structure-space Davidson. All quantum calculations and
builds ran through Hanhai22 Slurm, with 32 physical cores, cores/close,
one BLAS thread and 64 GiB requested. Historical executables were not rerun.
The unrelated local RI changes were neither transferred nor modified.

## Outcome: keep the previous fast production solver

**No new CI inverse model from this round is retained in production.**
The three experiments did not reduce the final 19-action CI solve.
The full-search-space variant markedly increased CI work. Their callbacks,
public API changes, auxiliary kernels and temporary unit tests were removed;
there is no production flag or fallback retaining them. Production numerical
source is again exactly revision 092dede3, apart from unrelated pre-existing
local RI work. The known 4-step/44-HVP/104.286136-s result remains the retained
baseline; it was not rerun after restoration.

What is retained is the algebra/evidence, a bounded diagonal-and-coupling
diagnostic, its offline reader, a deterministic arbitrary 600-structure
fixture, and the expanded existing 13-test regression launcher. Restored
production qualification 1265340 passes 13/13 tests in 6.02 s. The final
diagnostic compiles and runs successfully in 1265359 and 1265361.

## Mathematical experiments

Let R be the retained positive spin-energy inverse model and A a matrix of
constraint covectors. The exact inverse of the same positive model restricted
to ker(A^T) is R-R A(A^T R A)^-1 A^T R. This algebra is valid for arbitrary
independent selected columns and does not require orbital closure. The dense
restricted-model and nonorthogonal multiroot tests passed during qualification
of all three experiments. Valid algebra alone does not ensure a better
approximation to the physical shifted Hamiltonian inverse.

The tested choices were:

- A=M C: current requested-root constraints, using cached overlap images.
- A=[M C,(H-E M)v_next]: also use the nearest unrequested Ritz root already
  computed for the gap estimate. Its homogeneous correction equation follows
  from the projected Newton equation and Galerkin orthogonality.
- A=M V: the complete current Davidson metric basis. A bounded panel sweep
  formed A^T R A without retaining R A; a second vector inverse completed the
  constrained image. No extra physical H/M action or precision change was used.

The third model adds O(n_str*m^2) projection work for Davidson dimension m,
but also needs more Davidson directions in this case. It was not repaired by
changing a cutoff, suppressing a residual or adding a fallback.

A balanced coarse inverse using only the existing Davidson V is also ruled
out as a new-direction mechanism: V^T r=0 makes its difference from R r lie
in span(V), which the subsequent orthogonalization removes. This observation
does not invalidate coarse corrections in other linear solves.

Derivation and the explicitly unimplemented next carrier-inverse candidate
are recorded in article/source_column_hamiltonian.tex, Section 10.10. The PDF
was rebuilt and inspected. The carrier candidate is only a preconditioner,
not an exact inverse of the arbitrary selected physical metric; closure is
not asserted and no complete carrier allocation is proposed as a final design.

## Same-point cold CI and native results

The current-operator cold comparisons all use byte-identical endpoint input,
SHA256 9b226705643aca7915000e09e0afacf4642d0bf3ddf41809c80f1e17423417d3.
Seeds, exact actions and solve targets are unchanged. Only the inverse model
changes; this is not a repeated old-code SCF optimization. The 81-action
control is reused from the preceding diagnostic, not rerun.

| Inverse model | Cold CI actions | Active CI energy / Eh | CI residual / Eh | CI-only wall / s |
|---|---:|---:|---:|---:|
| Retained spin-energy, saved 1265254 | 81 | -18.662032474014442 | 9.7058e-8 | 7.0521 |
| Root constrained, 1265281 | 82 | -18.662032474014442 | 7.6887e-8 | 7.2092 |
| Nearest Ritz constrained, 1265311 | 82 | -18.662032474014456 | 8.2243e-8 | 7.6412 |
| Full projected inverse, 1265327 | 203 | -18.662032474014442 | 9.6424e-8 | 51.8722 |

These active CI energies are not total molecular energies. Cold counts must
not be substituted for the final warm solve's 19 actions.

| Native ION(0:2) | Steps | HVP | Final warm CI actions | Total energy / Eh | Application / s | Peak RSS / GiB |
|---|---:|---:|---:|---:|---:|---:|
| Retained saved 1265266_1 | 4 | 44 | 19 | -399.397232509245 | 104.286136 | 2.823959 |
| Root constrained 1265280_1 | 4 | 44 | 20 | -399.397232509244 | 99.978798 | 2.982868 |
| Nearest Ritz 1265309_1 | 4 | 44 | 21 | -399.397232509249 | 96.278069 | 3.175655 |
| Full projection 1265325_1 | 3 accepted | at least 41 | unavailable | not converged | cancelled | unavailable |

The first two trials show smaller single-run wall times, but no action-count
reduction and higher RSS. The retained and root-only runs use anode65,
affinity 96-127; the nearest-Ritz and full-projection runs use anode06,
also affinity 96-127. Both node types are AMD EPYC 9654. Do not attribute the
8% cross-node difference to an established algorithmic acceleration.
Native timing excludes profiler collection. GNU process wall times and
all precision/status fields are recorded separately in native_comparison.csv.
Cancellation leaves total HVP/time/RSS incomplete; accepted counters are
lower bounds, not completed-work statistics. No cancelled result is a speedup.

## Not an ION-specific test or implementation

The unchanged five archived decks are checked byte for byte before a new
fixture is added. The new nested 600-column arbitrary selection preserves
the first twenty archived columns and has 108 covalent, 161 singly ionic,
165 doubly ionic and 166 triply ionic Rumer structures. Its independently
variable structure coefficients are not enlarged to a closed/full space.

| Nearest-Ritz trial / 1265310 | Structures | Steps | HVP | Final energy / Eh | Application / s | Peak RSS / GiB |
|---|---:|---:|---:|---:|---:|---:|
| ION(0:1) | 5676 | 5 | 41 | -399.323427907218 | 42.214855 | 2.665840 |
| Arbitrary mixed20 | 20 | 4 | 38 | -399.095231633372 | 20.755742 | 2.609932 |
| Arbitrary mixed600 | 600 | 4 | 31 | -399.175606910787 | 21.103057 | 2.632000 |

These are converged trial results, not timings of a newly retained production
algorithm. There is no saved previous-version mixed600 timing, so no speedup
is claimed for it. The full-projection trial's corresponding results are also
retained in the CSV and evidence. Its mixed600 final CI count falls from six
to five, but time increases from 21.10 to 21.67 s: fewer actions alone is not
the end-to-end objective.

## What the exact-column audit actually shows

The retained production Jacobi inputs sum spin-product self terms. They are
not the complete structure diagonal, which includes interference terms.
The audit samples the first and last structure of each present ionic order;
one shared bounded H/M block action obtains exact columns. The model and exact
shifted diagonals are compared after division by G_II, matching the model's
coefficient normalization rather than comparing unrelated Rayleigh quotients.

For ION(0:2), the largest relative shifted-diagonal discrepancy in these six
samples is 28.69%, not an order-of-magnitude error. For the arbitrary mixed600
endpoint the largest of eight discrepancies is +91.37%. This is evidence of
nonuniform inverse-model quality, not an error in the exact forward action or
proof that diagonal error alone explains the CI iteration count.

Cross-occupation fractions use the norm sqrt(l^T R l), l=(H-E M)e_I,
with the retained occupation-block R. Complete occupation groups are mutually
orthogonal in this norm. It is a diagnostic model norm, not the physical
wave-function metric. Individual signed row products are never treated as
probabilities or thresholded. In the six ION columns, the cross-occupation
norm fraction ranges from 0.4163 to 0.7399. For the two doubly ionic columns,
the cross-ionic-order fraction is only 0.00531--0.00581. Thus most of their
cross-occupation coupling connects different arrangements within ionic order
two; merely merging scalar ionic-order blocks does not make a cheap inverse.

Occupation transfer distance is half the L1 change in spatial orbital
occupancies, not a determinant excitation rank. In the ION samples, distance
greater than two contributes at most 2.891% of the inverse norm. For the
doubly ionic sample patterns, distance<=2 covers 7146/47256 possible rows,
not the entire space. However the covalent ION patterns cover all rows at
distance<=2 simply because this input excludes higher ionic orders.

**The arbitrary-subset test prevents a false generalization:** one mixed600
covalent probe has a distance>2 tail of 21.13% of this norm, whereas other
probes have smaller tails. Consequently the ION observation is not a general
strict-sparsity theorem. No distance cutoff was added to the exact action.
The unimplemented carrier-inverse formula in the theory document retains
all occupation mixing instead of assuming that these tails vanish.

Both probe archives retain H, M, shifted columns, inverse images, occupation
bit patterns and sampled indices. tests/mcscf/ci_probe.py computes complete
distance-band statistics offline, with binary hashes, so a new statistic does
not require another electronic calculation. probe_analysis.json and
probe_mixed600_analysis.json are generated, not manually transcribed.

## Provenance

Qualification jobs 1265277, 1265300 and 1265318 each pass 13/13 tests for
their experimental numerical code; restored qualification 1265340 also
passes 13/13. Frozen bundle/source/build manifests, native outputs, CIAH
traces, GNU-time reports and probe payloads are in evidence.tar.gz. Exact
executables remain on Hanhai22. The saved controls are in the preceding
committed evidence directories and were not recomputed.

Preparation 1265285 failed because a test-data directory already existed,
not because of a quantum-chemistry failure. Its dependent array 1265286 was
cancelled without running. The generator was made idempotent; preparation
1265288 and subsequent bundles verify all original deck hashes. Only the
regressing native 1265325_1 was stopped during calculation. Accounting retains
all these states and missing fields. No task from this round remains running.
