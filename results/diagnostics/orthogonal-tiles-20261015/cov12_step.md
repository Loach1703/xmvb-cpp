# One physical SCF orbital step: CAS(12,12) complete COV

Hanhai22 job 1262608, anode20, production source 03919978. The test uses
`testdata/vbscf/csd_LOFLEA_vb.xmi`: original geometry and initial orbitals,
HAO, cc-pVDZ, libcint packed-exact integrals, 132 covalent structures.
Only ISCF is set to 9 (NEO) and the external macro-iteration limit to 1.
Davidson, analytic Hessian diagonal, gradient 1e-3 and energy 1e-7 Eh are
unchanged. No older executable was rerun.

| Quantity | Measured result |
|---|---:|
| Accepted outer steps | 1 |
| Initial energy / Eh | -398.533053417516 |
| Energy after the step / Eh | -399.074297731866 |
| Energy decrease / Eh | 0.541244314350 |
| Orbital HVP actions | 6 |
| Structure response actions | 6 |
| CIAH micro-iterations / keyframes | 10 / 3 |
| Rejected trials | 1 |
| Input preparation / s | 3.170278 |
| Program SCF-stage wall / s | 48.285492 |
| Program end-to-end wall / s | 51.464063 |
| GNU time process wall / s | 51.62 |
| GNU time process peak RSS / KiB | 2634968 |
| GNU time CPU use | 759% |
| Last CI Davidson iterations / H/S block actions | 99 / 99 |

This is not a convergence test: termination is `max_iterations`, and the
application exit code 2 is the expected one-step-limit status. The batch
preserves that status in `run.exit`; batch completion is not interpreted as
SCF convergence. The SCF timer includes initial objective evaluation and final
SCF diagnostics/export; it is not an exclusive timer for the accepted update.
The last CI residual is 7.6620e-8 Eh.

The process is pinned to CPUs 80--87, one LLC domain, with OMP_NUM_THREADS=8
and single-threaded BLAS. Sixteen CPUs are allocated so the helper can select
a complete eight-core cache domain. The application's `CPU threads: 16` field
reads the allocation rather than the OpenMP setting; actual affinity and
environment are retained. CPU utilization corresponds to about 7.59 cores.

The 0.506 s constructed CAS12 action/derivative fixture in the GitHub
comparison is not a timing model for this molecular SCF step. The molecular
step includes AO/core work, iterative CI solves, preconditioner/model work,
multiple HVPs, and keyframes. Counts alone do not attribute the 48.3 seconds
to any one component; no new hotspot attribution is claimed without stage
timing. This run does not compare molecular SCF speed with the remote source.

`cov12_step_evidence.tar.gz` retains the input, native stdout/stderr, trace,
GNU time record, build configuration, source/binary/input hashes, affinity,
node, and Slurm accounting. The runnable batch script is
`tests/mcscf/cov_step.sbatch`; the frozen binary remains in the cluster job
artifact.
