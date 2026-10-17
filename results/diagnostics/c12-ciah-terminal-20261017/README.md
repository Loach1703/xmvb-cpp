# ION(0:2): CIAH terminal-trial ordering

2026-10-17. Exact Libcint/cc-pVDZ, LOFLEA HAO CAS(12,12),
47,256 selected structures, ISCF=9/CIAH, structure-space Davidson.
Hanhai22/anode42, 32 physical cores, OpenMP cores/close, BLAS=1,
64-GiB Slurm allocation. No old executable was run for this change.

## Fix and unchanged accuracy

The old macro gate required a strictly positive computed energy reduction
before completing the trial gradient. A nearly stationary equal-energy or
slightly higher-energy candidate was therefore rejected before the ordinary
two-condition stopping rule could inspect it. The preceding ION(0:2) record
exhibited this ordering: source -399.39723250925095 Eh, last rejected trial
-399.39723250925084 Eh, after five accepted steps with a small true gradient.

The new ordering is:

1. Normal downhill candidates retain their existing acceptance path.
2. For a non-downhill candidate whose absolute energy change is below the
   **unchanged 1e-7 Eh** threshold, complete its cached forward evaluation
   to obtain the **true relaxed gradient**.
3. Apply the original gradient predicate, at that candidate's own orbital
   point, with the **unchanged 1e-3** threshold. This includes the existing
   physical dual-norm requirement for complete OEO, when applicable.
4. Accept this candidate as terminal only if both original conditions pass.
   Otherwise retain strict descent/backtracking. A small model gradient or
   the source gradient cannot substitute for the candidate's true gradient.

No numerical-energy slack, new tolerance, gradient-only stopping, stability
certification, CI approximation or physical-space change is introduced.
The new check is not performed on ordinary downhill steps. It neither adds
an HVP nor rebuilds integrals. Its completed gradient, geometry and stopping
verdict are retained, so the outer loop does not repeat their construction
or certification. The redundant `AcceptedCiahKeyframe` wrapper was removed;
the driver now consumes `CiahMacroStep` directly. Accepted and trial points
share one geometry builder and one gradient predicate, without mutating
the accepted objective or copying its molecular integrals.

## Results from native runs

All three input files are byte-identical:
`48b44d64b773d4abaf388fcf728d478a09ab76f5426f889822b1dbc8e031432f`.
The original initial orbitals, ITMAX=2000, numerical thresholds and memory
budgets are unchanged. Both controls below are saved results, not reruns.

| Version / job | Status | Accepted steps | Total HVP | Final energy / Eh | Application total / s | Peak RSS / GiB |
|---|---|---:|---:|---:|---:|---:|
| Historical / 1265093_1 | Converged | 6 | 86 | -399.397232509251 | 1631.498160 | 3.583797 |
| Spin-block preconditioner, before fix / 1265114_1 | Terminal rejection | 5 | 86 | -399.397232509251 | 173.676211 | 3.066448 |
| Current terminal ordering / 1265186_1 | **Converged, exit 0** | **6** | **86** | **-399.397232509252** | **153.175347** | **2.886070** |

The failed control's total includes 83 accepted-step HVPs and three in the
failed sixth macro. The current per-step counts are **12,12,14,18,27,3**.
There are no rejected trials in its accepted trace. Final projected gradient
infinity norm is **4.00347653e-7**, and the last energy change is
**-1.81898940e-12 Eh**. Both original stopping conditions pass.

The current energy differs from the historical result by about -9.66e-13 Eh.
The observed total-time ratio to that historical record is 10.65x; the major
speedup is the preceding spin-block preconditioner, not a claim that this
termination repair alone accelerates the physical kernel by that factor.
The preceding fast failed record took 173.68 s versus 153.18 s now; neither
its failure nor its time should be mislabeled as a converged control. These
are single measurements on the same node, not repetitions on an exclusive
allocation with an identical CPU set.

Time in the table is the application's end-to-end timer. GNU-time complete
process times are 1634.85, 176.46 and 156.06 s respectively, including trailing
output/teardown. RSS is GNU-time process peak, not sampled Slurm MaxRSS or
the allocated memory. Current CPU affinity is 96-127 (32 cores). Slurm job
and native step both completed with exit 0.

## Verification scope

Qualification **1265176** passed **11/11 tests in 3.64 s**, including new
terminal-candidate tests in `test_keyframe`. These cover equal energy,
one-ULP increases/decreases, an energy change beyond tolerance, the strict
energy-threshold boundary, a nonconverged true gradient, and a missing true
gradient, in both complete and selected chart modes. The remaining tests
cover CIAH, selected metric/transport, coupled derivatives and direct-CI.

The current native run's last computed energy change was negative. Thus it
demonstrates complete convergence with the new implementation, not a
bitwise reproduction of the previous opposite-ULP rejection. The directed
unit tests explicitly exercise nondecreasing terminal candidates; no old
calculation was repeated merely to reproduce a rounding sign.

The frozen executable was produced by qualification 1265176 and copied by
**1265185**; its SHA256 is
`6b88660f308f303d22746f5c2e7d217edba77f1f61761475ad53a08b4a12a75f`.
Its source manifest and Release CMake configuration are in the input bundle.
The execution mirror excludes the unrelated local RI edits. Native run
**1265186_1** uses the frozen executable, not a subsequently rebuilt binary.

## Saved records

- `comparison.csv`: programmatically parsed three-version comparison, including
  input hashes, accepted/failed HVP counts, last energy changes, true projected
  gradients, process/application times and RSS.
- `accounting.txt`: Slurm qualification, freeze and native-job accounting.
- `evidence.tar.gz`: frozen input/build metadata (without executable) and the
  complete current native output, trace, orbitals and GNU-time record.
- Historical records: `c12-ionic-metric-20261017/evidence.tar.gz` contains both
  earlier ION(0:2) runs. No controls are newly executed or duplicated here.

The live current directory is
`/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015/c12-spaces-1265186/ion02`
on Hanhai22. Existing `tests/mcscf/c12_spaces.py` parses all three records;
select its `ion02` row for the table. The three roots are `run-1265093`,
`run-1265114` and `run-1265186` after extracting the two evidence archives.
