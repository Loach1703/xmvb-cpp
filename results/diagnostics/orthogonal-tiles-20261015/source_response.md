# Bounded selected-structure coefficient responses

## Change and mathematical contract

`StructureAction::apply_response()` now compares two exact contraction orders
for its complete H/M image and adjoint request. The factor order expands the
response only on original touched spin products and reuses those values in
the shared source projection and integral/residual-generator adjoints. It
does not allocate an orthogonal response, sigma, projected-sigma copy, or
response residual product array. The original support is expanded once.

For Y = T W C and Q = T W q, the integral contribution is 2 J(Q,Y). The
generator contribution is **both** 2 D((H-E)Q,Y) and 2 D((H-E)Y,Q). Neither
residual is projected into the selected structure space before contraction.
State weights and the separate off-diagonal normalization/orthogonality
multipliers retain the existing gauge convention. Section8.4 of the TeX/PDF
derives the equations and leading operation counts.

The arithmetic comparison includes one integral adjoint, two residual
generator adjoints and the additional H/M image sweep. It has no fitted
constant or accuracy threshold, and never changes the selected space.
Broad-support cases can still use the exact product contraction order.

Ownership remains small and explicit: `source_response.cpp` owns the new
consumer; `source_tiles.cpp` shares its support-column projection with forward
H/M; `mcscf/ci` owns the numerical factor kernels. The new consumer does not
implement a second spin scatter, Hamiltonian or exterior-transform kernel.

## Verification and execution

All builds and workloads ran in Hanhai22 Slurm compute allocations, anode26,
8 bound cores, GCC14.1.0/CMake3.30.5, existing Release O3 flags and BLAS1.
Jobs1261034 and1261052 each passed all9 CTests. Production numerical code is
the same across those jobs;1261052 adds the CAS14 benchmark fixture/CLI only.
Inputs, flags, source/binary hashes and raw timing logs are retained.

New small-system tests cover positive/negative/no common exchange parity,
unequal alpha/beta strings, unequal state weights, off-diagonal gauge terms
and exactly zero response. They clear accepted orthogonal/residual arrays
before calling the new consumer, proving that it does not secretly borrow
those arrays. H/M is compared with independently expanded complete-carrier
sigma; the adjoint is compared with the separate directional-gradient path.

The larger kernel fixtures use deterministic synthetic active integrals,
not complete molecular optimization. CAS14 uses **two different covalent
structures** and projects Q against C in the full structure metric, then
normalizes its nonzero physical tangent. CAS12 single-structure response is
an algebra/plumbing benchmark, not a nontrivial horizontal CI tangent.

## Native response measurements

Times below cover `apply_response`: both structure H/M images and the
orthogonal integral/generator adjoints. They exclude accepted-state preparation
and the final small active-integral pullback. RSS is the native process peak,
including accepted preparation/caches and both cold/warm calls. It is not
the isolated response arena size. Native runs exclude the reference oracle.

| Active space | VB structures | Exact order | Cold s | Warm s | Native RSS MiB |
|---|---:|---|---:|---:|---:|
|12/12|1|source|0.09307|0.07000|87.45|
|12/12|3|product|0.21499|0.23242|100.73|
|12/12|20|product|0.22573|0.21617|100.46|
|12/12|132|product|0.22152|0.20987|100.88|
|14/14|2|source|2.85399|2.83700|293.02|

The CAS12/CAS14 source-image arenas are bounded by63,866,880/52,880,256 bytes,
both below the shared64 MiB factor-workspace ceiling. Original-support input
and output vectors, small returned adjoints and immutable metadata are
separate. Every factor numerical workspace is shared by the OpenMP team;
there is no per-worker product array.

An initial new-code variant (job1261024) inherited the larger general H/M
workspace: CAS12 source-image bound113,541,120 bytes, warm0.12343 s and
101.33 MiB RSS. It was corrected before the final version by sharing the
existing64 MiB factor ceiling, not increasing a budget to obtain speed.
Its observations are retained but are not historical-baseline runs.

Separate reference processes retain the current full-product contraction,
not an old executable. Their **adjoint-only** call took0.18863 s for CAS12
one structure and4.00804 s for CAS14 two structures. That reference includes
the final pullback but no H/M images: the timing scopes differ, so no exact
whole-consumer speedup ratio is claimed. The source consumer already
includes both images and all residual terms.

CAS12 source H-image/adjoint errors are4.44e-16/1.49e-15; CAS14 errors are
4.00e-16/1.51e-15. Product-order cases retain their independent H-image
checks and agree with the separate response adjoint. The CAS14 native process
used8.83 s in total,768% CPU and293.02 MiB RSS. Its reference process peaked
at455.14 MiB, but also includes extra independent oracle allocations; that
is not a matched production-memory baseline.

## Complete molecular regression versus saved data

Exact libcint, NEO/Davidson, analytic Hessian diagonal, projected-gradient
threshold1e-3 and energy-change threshold1e-7 Eh are unchanged. Saved baseline
is job1260422 on the same node/eight cores; no historical executable was rerun.
Final current application job is1261034.

| 240 space | Steps | HVP | Saved/current SCF s | Saved/current RSS MiB | Current energy Eh |
|---|---:|---:|---:|---:|---:|
|COV|3|20|6.676/6.447|1528.6/1509.4|-343.239211501262|
|ION(0:1)|3|27|8.761/8.697|1517.8/1547.4|-343.407736636063|
|ION(0:2)|4|32|10.834/10.950|1518.7/1514.7|-343.443415326850|
|ION(0:3)|4|32|10.983/11.089|1515.6/1499.4|-343.446256564846|
|FULL|4|30|9.623/9.749|1525.7/1531.2|-343.446302545203|

Steps/HVP are unchanged, maximum final-energy difference2.96e-12 Eh.
Whole-task time/RSS remains essentially unchanged; these broad-support
response calls do not use the new source order. No end-to-end molecular
speedup is inferred from the smaller-source consumer measurements.

## Remaining work

This is one complete consumer migration, not a bounded-memory complete SCF.
Accepted-state preparation still retains full orthogonal/residual arrays.
Orbital-direction actions, directional generator adjoints, the selected-state
producer and restricted metric metadata still require migration. The new
consumer calls its integral and two residual adjoints separately; fusing
their shared moments remains possible and is not claimed implemented here.

Data: `source_response.csv`, `source_response_240.csv`, `source_response.sacct`,
raw `source-response-1261034`, `source-response-1261052` and the final molecular
`source-response-scf-1261034` records. Parsing is reproducible with
`tests/mcscf/summarize_factor_gradient.py`; no manual extraction of trajectories.
