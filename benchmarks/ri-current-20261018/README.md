# Current RI native measurements

## Scope and reproducibility

All five cases converge with CIAH/Davidson on Hanhai25, exact same original
input files and orbital guesses as their saved exact references. Only the
explicit CLI integral mode changes to `--standard-two-electron-mode ri`.
Primary/auxiliary bases are cc-pVDZ/cc-pVDZ-JKFIT for241,C12,C14 and
def2-TZVP/def2-TZVP-JKFIT for FeCl2,MnF2. Outer thresholds and the4GiB shared
CI allowance are unchanged. RI energies are not required to equal exact
energies: their signed differences are reported, not treated as zero.

The frozen binary is the current qualified3c31ec4c source, from qualification
293572, SHA256
`cfc05fc70ff58d98ce03066e8e3112202224b5fa6ad09128aeb191a2b3ddb390`.
No rebuild or production-code change was needed. The pre-existing uncommitted
changes to `ri_operator.cpp` and its test are **not** part of this frozen binary
and remain untouched. This benchmark must not be attributed to that patch.

Native job293579 retains one32-core allocation across separately measured
processes, preventing overlap between measured cases. Node anode055,
CPU-256C768GB, AMD EPYC9755, CPUs128-159; OpenMP32, cores/close, BLAS1.
All build/workload execution follows the Slurm compute-node policy. GNU
`time -v` measures complete-process wall time and maximum RSS; mean cores are
(user+system)/wall, which includes OpenMP spinning/waiting.

Saved exact controls are on other nodes of the same cluster/partition.
241 and C14 use the near-preceding5b46ef5d source, C12 uses3c31ec4c.
FeCl2 and MnF2 are older retained-policy references from291714, not matched
current-source controls. Input hashes match all five. Therefore the table is
a comparison against saved observations, not a same-node isolated integral
backend speedup experiment. No historical executable was rerun.

## Converged native observations

| Case | Exact steps/HVP | RI steps/HVP | Exact / RI wall s | Exact / RI RSS GiB | RI mean cores | RI minus exact mEh |
|---|---:|---:|---:|---:|---:|---:|
|241|3/20|3/20|1.26 /0.99|0.971 /0.189|29.04|+0.309736|
|C12 COV132|4/33|4/33|12.07 /12.26|3.435 /1.427|27.04|+0.370154|
|C14 COV429|4/46|4/46|168.88 /165.24|8.139 /4.033|25.95|+0.555423|
|FeCl2 OEO|13/170|12/229|3.68 /7.67|0.434 /0.245|31.08|+0.096803|
|MnF2 HAO|11/160|9/132|2.66 /3.56|0.495 /0.185|31.23|+0.079013|

Memory reductions are substantial, especially for the organic cases; speed
is not uniformly better. C12 and C14 retain the same per-macro HVP counts
6:16:8:3 and5:10:28:3 as exact. FeCl2 uses more HVPs despite one fewer macro;
MnF2 uses fewer of both, but is still slower than its older saved exact run.
Full RI energies, input/binary hashes, action counts, timings and reference
provenance are in `results.csv` and `exact-references.csv`.

## C12 attribution, without changing the integral implementation

Separate phase job293581 uses the same binary,node and core mask as native
293579. It also converges in4 steps/33 HVPs. Its times are diagnostic scopes,
not a replacement for native measurements. `c12-phases.csv` compares them
with saved exact audit293573 on anode003, also32 cores and the same partition.
Nested scopes must not be summed.

| Scope | Saved exact s | RI s | Calls exact / RI |
|---|---:|---:|---:|
|HVP apply|3.412447|4.533356|33 /33|
|AO response|0.439302|1.282917|33 /33|
|Active integral direction|0.005150|0.000531|33 /33|
|Two-electron response|0.082991|0.160238|33 /33|
|CI sigma|2.604999|2.510935|79 /79|
|CI response|2.343669|2.226024|47 /47|

The RI loading path avoids building raw exact AO ERIs and their pair graph.
The orthogonal direct-CI operator still consumes the same active pair kernel,
so changing the AO integral backend does not reduce spin-product contraction
work. These observations locate the additional response cost on the AO side,
not in the selected-space sigma scheduler.

Source review identifies concrete follow-up candidates, not measured individual
function costs: the first adaptive AO-H1E action at a point runs both dense and
spectral implementations for calibration. The spectral route copies its dense
source/adjoint and decomposes them separately before two RI sweeps. The code
already has an explicit-factor overload, but this HVP caller supplies matrices.
The next derivation should identify the exact low-rank source/adjoint factors
already available upstream, then fuse their shared accepted factors and scans.
This is not a reason to relax spectral cutoffs or change auxiliary bases.
Likewise, independently constructed RI one-index response contexts may share
accepted data only after verifying their orbital point and lifetime.

## Files and replay

Raw artifacts are under
`results/diagnostics/ri-current-20261018/{native-293579,phase-293581}`;
the remote frozen binaries, complete outputs and Slurm accounting remain on
Hanhai25. `run.sbatch` stages the original inputs and fitting bases, records
hashes/affinity/resources, and supports native or separate phase mode.

Recreate the CSV from saved logs, without calculating:

```bash
python3 benchmarks/ri-current-20261018/summarize.py \
  results/diagnostics/ri-current-20261018/native-293579 \
  benchmarks/ri-current-20261018/exact-references.csv \
  benchmarks/ri-current-20261018/results.csv
```
