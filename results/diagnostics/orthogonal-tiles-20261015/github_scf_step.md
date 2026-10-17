# Molecular COV/FULL one-step comparison with GitHub

## Setup

- Baseline: unmodified GitHub source 946934570c096ba5072923effa7fb063736731f1.
  Its production executable was built in Hanhai22 Slurm job 1262620, reusing
  the already-qualified Release library build; no production code was patched.
- Current: frozen executable from job 1262608, numerical source 03919978
  (e5ef2fe1 changes only the test runner/report).
- Original LOFLEA geometry, orbital guess and HAO supports, cc-pVDZ,
  libcint packed-exact integrals, 12 electrons in 12 active orbitals.
  Set STR=COV or FULL and ISCF=9; NEO/Davidson, analytic Hessian diagonal,
  gradient 1e-3 and energy 1e-7 Eh. Only one outer accepted step requested.
- Jobs 1262622 (COV) and 1262623 (FULL) ran concurrently on anode50. Within
  each job the variants ran sequentially on identical cores: COV CPUs 80--87,
  FULL CPUs 176--183, one eight-core LLC domain each. OpenMP 8, BLAS 1;
  16 CPUs allocated per job to admit a complete LLC domain, 48 GiB memory.
  These are single-run measurements on a shared node, not repeat averages.
- Each variant/space was executed once in this matched campaign. The earlier
  current-only COV job 1262608 is not substituted into the matched table.

## Results

| Space | Version | Structures | Accepted steps | Orbital HVP / response | SCF stage (s) | Native process wall (s) | Peak RSS (GiB) |
|---|---|---:|---:|---:|---:|---:|---:|
| COV | GitHub | 132 | 1 | 6 / 6 | 53.464025 | 58.06 | 2.61 |
| COV | Current | 132 | 1 | 6 / 6 | 48.055195 | 52.66 | 2.52 |
| FULL | GitHub | 226512 | 1 | 8 / 8 | 74.965957 | 93.77 | 3.37 |
| FULL | Current | 226512 | 1 | 8 / 8 | 112.829673 | 132.88 | 3.54 |

| Space | Initial energy, both (Eh) | GitHub after step (Eh) | Current after step (Eh) | Absolute endpoint difference (Eh) |
|---|---:|---:|---:|---:|
| COV | -398.533053417516 | -399.074297757091 | -399.074297760467 | 3.38e-9 |
| FULL | -398.880685140388 | -399.322429054264 | -399.322429045702 | 8.56e-9 |

All four application runs terminated at the requested `max_iterations=1`,
with exit code 2; none is labeled converged. The final CI solve uses 99 H/S
block actions for COV and 39 for FULL, in both versions. Both versions use
10 micro-iterations / 3 keyframes in COV and 13 / 3 in FULL.

COV SCF-stage time decreases 10.1% and peak RSS decreases 3.6%. FULL SCF-stage
time instead increases 50.5%, and RSS increases 4.8%. The reported improvement
of the constructed fixed-point kernel therefore does not establish a general
molecular SCF speedup. FULL is a measured regression and remains unresolved.
The identical action counts rule out a larger recorded HVP/CI-action count
as its explanation. They do not locate the expensive function: per-action
work, data movement, synchronization, and other SCF work require attribution.
Native CPU time for FULL increases from 446.71 to 763.49 CPU-seconds; this
includes OpenMP busy-wait time and is not a FLOP count.

## Timing scope and evidence

The application's SCF-stage timer includes initial objective evaluation and
final SCF diagnostics, not just the accepted orbital update. The native GNU
time includes the complete process and output writing. In particular, each
FULL stdout is about 158 MiB because it prints all structure weights; program
end-to-end timers (80.97 / 119.05 s) do not include all of that final output.
Native process wall is therefore also retained rather than silently replacing
it with the application's shorter timer. RSS is GNU time process peak.

`github_scf_step.csv` is parsed from saved outputs by
`tests/mcscf/summarize_scf_step.py`, without launching workloads. Inputs, COV
stdout, FULL diagnostic extracts, trace, GNU time, binary/input hashes,
affinity, configuration and accounting are in
`github_scf_step_evidence.tar.gz`. The two unabridged 158 MiB FULL stdout
files remain in the local and cluster artifacts, with SHA-256 hashes in the
archive, rather than adding hundreds of thousands of printed weights to git.
Exact executables remain in the cluster job artifacts. No optimizer or
numerical kernel was changed for these tests.
