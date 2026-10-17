# Independent TNVB: build, run and numerical acceptance

## Master synchronization: 2026-10-17

Master `03d469cc` is integrated with TNVB restricted to ISCF=8/P-LBFGS.
The C++23 build includes the production CLI and all TNVB numerical tools.
Five focused CTests pass: input invariants, gradient adapter, RDM policy,
acceptance-summary logic, and trial isolation with an accepted orbital step.
The normal benzene ENERGY CLI run and fixed-orbital exact/RI energy-gradient
comparisons also pass; maximum same-orbital gradient differences are
5.11e-15 / 4.94e-15. This is not a new full C14 optimization acceptance.
The isolated verification build is
`/pool1/home/chentz/apps/xmvb4.0/build/cpp-tnvb-master-20261017`.
Its LAPACK and LAPACKE explicitly use the previously verified dynamic
OpenBLAS library recorded in the older `cpp-tnvb/CMakeCache.txt`; automatic
selection of the older static library caused SIGILL on this host.
No existing installation or submission script is replaced.

## Scope

The `tnvb` branch supplies a separate `xmvb_tnvb` library. Production uses
singlet VB structures, SU(2) Block2, fixed bond dimensions, DIRECT
total-state construction and FULL 3-RDM. HAO and OEO use the same electronic
adapter. Selected structures remain a selected space, not a complete CAS.
The default ordinary determinant path and its second-order caches remain
available; a TNVB evaluation creates neither determinants nor an
`AcceptedPointContext`. TNVB uses ISCF=8, the shared nonredundant orbital-block
P-LBFGS optimizer (`Plbfgs`), after merging master `03d469cc` on 2026-10-17.
It uses the same first-order electronic adapter and does not require an
`AcceptedPointContext`, a determinant HVP, or a new input keyword. Its original
stopping metric is the projected reduced-gradient infinity norm. Enabling TNVB
does not change that optimizer's numerical stopping thresholds or algorithm.
Master removed native ISCF=5; this branch does not restore or alias it.

Only these new input keywords are registered:

| Keyword | Rule |
|---|---|
| `TNVB` | Enable the independent electronic backend. |
| `TNBOND` | Required positive integer; cap for each structure's SU(2) MPS. |
| `TNSUMBOND` | Required for optimization, positive and at least `TNBOND`. |

For ENERGY, use `VBSCF ISCF=8 ITMAX=0 TNVB TNBOND=128`. No total MPS or
1/2/3-RDM is constructed. For optimization, use
`VBSCF ISCF=8 ITMAX=200 TNVB TNBOND=128 TNSUMBOND=256`.
The CLI override is `--optimizer-backend plbfgs`. In ENERGY (`ITMAX=0`),
the fixed orbitals are evaluated without optimization.
These are additions to the existing `$CTRL`; geometry, orbital guess and
VB structures are still needed. Structures can be supplied explicitly through
`$STR` or generated using the existing supported `STR=...` structure classes,
including `STR=COV`. Both sources feed the same TNVB MPS construction without
determinant expansion. Final command-line `--max-iterations`
overrides `ITMAX`. Omitted `ITMAX` retains the project default 2000; specify
zero explicitly for ENERGY. TNVB defaults to dense H/S solving. An explicit
unsupported eigensolver or optimizer request is a capability error.

No old TNVB keywords, adaptive scheduling, contracted backend, TNHVP, CIAH,
state averaging, open-shell treatment or nuclear derivatives are registered
by this module. Generic unknown-keyword handling remains the existing
project behavior, without compatibility or old-keyword rejection branches.

## Build isolation

The project requires CMake >=3.20 and a C++23 compiler.
`XMVB_ENABLE_TNVB` defaults to OFF. This build does not require Python or
Block2. Configure an independent build directory and dependency roots:

```sh
cmake -S /path/to/xmvb-cpp -B /path/to/build-off \
  -DCMAKE_BUILD_TYPE=Release -DXMVB_ENABLE_TNVB=OFF \
  -DLAPACK_ROOT_DIR=/path/to/openblas \
  -DCINT_ROOT_DIR=/path/to/libcint -DEIGEN3_ROOT_DIR=/path/to/eigen
cmake --build /path/to/build-off --parallel 4
ctest --test-dir /path/to/build-off --output-on-failure
```

For ON, also configure matching paths once:

```sh
cmake -S /path/to/xmvb-cpp -B /path/to/build-tnvb \
  -DCMAKE_BUILD_TYPE=Release -DXMVB_ENABLE_TNVB=ON \
  -DLAPACK_ROOT_DIR=/path/to/openblas \
  -DCINT_ROOT_DIR=/path/to/libcint -DEIGEN3_ROOT_DIR=/path/to/eigen \
  -DXMVB_TNVB_PYTHON=/path/to/python \
  -DXMVB_TNVB_BLOCK2_MODULE=/path/to/compiled-block2 \
  -DXMVB_TNVB_BLOCK2_SOURCE=/path/to/block2-source \
  -DXMVB_TNVB_LIBRARY_DIR=/path/to/matching-gcc/lib64
cmake --build /path/to/build-tnvb --parallel 4
ctest --test-dir /path/to/build-tnvb --output-on-failure
cmake --install /path/to/build-tnvb --prefix /path/to/new-install
```

Use a BLAS built for the execution node, preferably dynamic CPU dispatch.
The original fixed SKYLAKEX OpenBLAS caused SIGILL on this host's AVX2 CPU;
the existing dynamic-architecture LP64 OpenBLAS solved that issue. GCC 13
and matched dependencies avoid the earlier incompatible-binary problem.
The local validated cache and exact dependency paths are under
`/pool1/home/chentz/apps/xmvb4.0/build/cpp-tnvb/CMakeCache.txt`.

Keep the matching compiler runtime active when building and running the
original project's test binaries. On this machine that is:

```sh
source /etc/profile.d/modules.sh
module load gcc/13.2.0
```

Without it, the shared ordinary `cli_options` test can load the shell's
Anaconda C++ runtime and segfault in both OFF and ON builds. The installed
TNVB executable and new acceptance tools have explicit matching runtime
paths; this shell requirement is not another TNVB input or worker setting.

Installation contains `bin/xmvb-cpp.exe`, `lib/libxmvb_tnvb.a`, its own
`libexec/xmvb-cpp/tnvb/{evaluate.py,pair_core.py}`, and basis assets in
`share/xmvb-cpp/basis`. Configured external Python/Block2/runtime dependencies
must remain accessible. This is not a self-contained binary distribution.
Existing project/local basis files retain precedence over installed assets.
Only the ON executable gains this independent install rule.

The build stages an owned worker copy; production never imports the legacy
XMVB worker. Installation resolves its worker relative to the executable.
Child-only Python paths and matched C++/OpenMP preloads do not change the
parent's environment; inherited preloads are retained. Each evaluator has
unique scratch/cache ownership and each request has its own files. A warm
state is only an initial guess: H/S and RDM are recomputed for current
integrals. Trial numerical results cannot replace the accepted point until
the optimizer explicitly commits them.

## Run and resources

```sh
/path/to/new-install/bin/xmvb-cpp.exe calculation.xmi > calculation.xmo 2>&1
```

RI is selected using the existing `RI`/`AUXBASIS` input conventions or
`--standard-two-electron-mode ri`. The worker receives active-space ERIs;
the C++ backpropagation retains RI factors and does not reconstruct AO ERIs.
RI and exact calculations are separate integral models: agreement between
the two electronic backends is tested within each model.

The example `src/tnvb/templates/submit.slurm` needs only the input and installed
executable. It runs in the submission directory. Adjust the Slurm CPU/memory
request for the system; no Python/worker/backend switches are needed.
Worker processes and Block2 threads are chosen from Slurm limits/CPU affinity,
the useful work and available Slurm memory. The actual state/H/S/RDM layout is
printed. `XMVB_BLOCK2_WORKERS`, `XMVB_BLOCK2_THREADS`,
`XMVB_BLOCK2_RDM_THREADS`, `XMVB_BLOCK2_STACK_GB`,
`XMVB_BLOCK2_RDM_STACK_GB` and `XMVB_BLOCK2_SCRATCH_PARENT` remain advanced
resource overrides, not ordinary input keywords. A Block2 stack is not the
whole job's memory budget. Multi-process copies and C++ AO/response work also
consume memory.

Finite `fit_dE>1e-7 Eh` emits a warning and continues. It is not a user-configured
convergence limit. Nonfinite data, invalid structure metric, worker failure,
eigenpair norm/residual failure and RDM check residuals above `1e-7` are errors.
Finite RDM check residuals above `1e-8` and at most `1e-7` emit a warning and
continue. Both limits are fixed in the C++/worker production path; inherited
`XMVB_BLOCK2_RDM_CHECK_TOL` cannot override them. This changes numerical
consistency handling, not orbital convergence parameters or RDM algorithms.
At finite bond dimensions, the gradient does not differentiate all truncation
operations. Check bond convergence before interpreting an approximate optimum.
ENERGY output reports no orbital convergence or gradient norms. Optimization
retains H/S, all existing structure weights, orbital and density analyses.
The CLI regression also checks Molden export in ENERGY and gradient modes,
all 24 orbitals of the benzene fixture, natural-orbital and population sections
after a full-gradient evaluation, and their absence in ENERGY output.
ISCF=8 integration checks and their explicitly bounded C14 scope are recorded
separately in `tnvb_iscf8_integration.md`; the original ISCF=5 acceptance archive
is preserved as historical evidence. The P-LBFGS optimizer's Armijo trials evaluate full
gradients, including FULL 3-RDM, even for rejected trials. Joining this optimizer
does not by itself reduce the cost of a gradient evaluation.

## Reproducible acceptance tools

Build optional numerical tools explicitly:

```sh
cmake --build /path/to/build-tnvb --target \
  compare_tnvb check_tnvb_directions optimize_tnvb \
  prepare_tnvb_oeo export_tnvb_request --parallel 4
```

All results are JSON at 17-digit precision, with explicit pass/fail, not tests
of rounded `.xmo` energies:

```sh
compare_tnvb INPUT REPORT.json exact gradient
compare_tnvb INPUT REPORT.json ri energy
check_tnvb_directions INPUT FD.json exact
optimize_tnvb INPUT OPT.json exact 200
```

`compare_tnvb` accepts an optional final-orbital plain-text table. It expands
determinants only in its explicit small-system reference branch, normalizes
both structure bases, aligns the root's overall phase, compares H/S,
coefficients and four weight definitions, and applies the common orbital
pullback for gradients. `check_tnvb_directions` checks independent S, h,
packed-g, inactive-orbital and active-orbital directions at three step sizes;
two adjacent passing step sizes are required.

`optimize_tnvb` uses the same guess, ISCF=8 and convergence settings for both
backends, retains failed/unconverged outcomes, compares both endpoints and
cross-evaluates their exported physical orbitals. A different stationary point
is not labeled as a passing optimization match. C14 is capped at 200 accepted
steps. Parent peak RSS and largest reaped child peak RSS are recorded separately;
neither is claimed to be the aggregate concurrent job peak. Slurm/cgroup
accounting is needed for that quantity.

The benzene profile retains its strict H/S, energy, coefficient, weight and
gradient gates. The user-confirmed C14 profile uses energy <=1e-6 Eh,
gradient <=1e-7, and coefficient/weight errors <=1e-4. Its finite-D H/S
differences are reported separately, including whether they satisfy benzene's
stricter matrix limits. Nonfinite matrices remain errors for both profiles;
The historical numerical acceptance profile retains RDM identity <=1e-8;
production now warns above 1e-8 and stops above 1e-7. Warning-band continuation
alone does not establish that stricter acceptance criterion. Eigenpair/metric
checks are never disabled. C14 final
optimization additionally requires both backends to converge under the
original optimizer criteria, fit <=1e-7 and exported-energy consistency.

`optimize_tnvb` uses only the original optimizer convergence settings:
energy tolerance 1e-7 Eh and raw gradient L2 tolerance 1e-3. There are no
additional tolerance arguments. The user explicitly requires unchanged
convergence parameters. Earlier reports labeled `stricter_diagnostic` are
retained solely as historical diagnostics and cannot supply formal acceptance.
C14 retains its original tolerances and 200-step limit.

The user-confirmed cross-backend optimized terminal-energy difference target
is 1e-6 Eh for both benzene and C14. This comparison target does not change the
optimizer's 1e-7 energy-change or 1e-3 gradient L2 convergence settings. Benzene
same-orbital energy/H/S/coefficient/weight/gradient gates remain unchanged.

`prepare_tnvb_oeo HAO.xmi OEO.xmi` exports the actual supported, normalized
HAO physical guess into full-AO OEO coefficients. The original autoVB GUS
contains repeated full-AO columns which are valid only after HAO support
projection; replacing the `ORBTYP` token alone gives a singular OEO guess.
The checked-in benzene OEO fixture includes this conversion.

`compare_legacy_worker.py` is test-only: explicitly provide a frozen legacy
worker and a request from `export_tnvb_request`. It compares H/S, root,
coefficients and all transferred density/contraction buffers. It has no
production import or fallback role.

`diagnose_overlap.py REQUEST WORKER REPORT` is a small-system (K<=14),
test-only batched-LU oracle for structure overlaps. Optional `--bond` and
`--pair-order right-endpoint` diagnose finite-D and pair-product-order effects
without changing the production input. Its explicit CSF enumeration belongs
only to the test reference; production TNVB remains unexpanded.

See `tnvb_gradient_conventions.md` for the mathematical interface and
`tnvb_acceptance_report.md` for the completed acceptance results and limitations.

`tests/tnvb/summarize_acceptance.py BUILD_EVIDENCE C14_EVIDENCE --output SUMMARY.json`
independently rechecks the saved numerical errors, finite-difference interval,
original optimizer criteria, endpoint gradients, C14 fit/export checks and
iteration cap. It never edits its input evidence or jobs. Missing reports are
`pending_evidence`, not proof that a job is running or successful. Original
optimization checks and supplementary stricter diagnostics are separate;
passing supplementary checks cannot turn an original failure into a pass.
Historical C14 matrix-only failures retain their original flag and measured
H/S errors while applying the explicitly user-confirmed C14 profile. The
old benzene optimization flags are likewise retained when reassessing their
measured errors against the user-updated 1e-6 terminal-energy target; changed
convergence settings and unexplained failures cannot be reclassified. The
summary is numerical evidence only, not a substitute for the final source,
build, installation, isolation and requirement-by-requirement audit.
