# Native C12/C14 follow-up

The user switched performance validation from C16 to C12/C14. New C16
jobs293254/293255 were cancelled before completion; their partial timing is
not an end-to-end result. Historical executables are not rerun.

Frozen executable: qualification293253,6/6 tests in6.57s, SHA256
`141283957128192fe63fac1a82120b1cb98202876b06242bf97bad5f2a81a7f5`.
It includes compact beta target rows and fused exchange weighting.

Use original archived inputs, exact Libcint/cc-pVDZ, CIAH/Davidson,32 physical
cores, OpenMP close/cores, BLAS1. All cases run to their existing convergence
criteria and original input iteration limits. No extra SCF step cap is set.
C12 here is the established CAS(12,12) LOFLEA benchmark, not a carbon-C12
molecule. C14 is anthracene CAS(14,14),COV429.

Cases are C12 COV132,ION(0:2)47256,FULL226512 and C14 COV429. The main runs
use the current4GiB shared CI budget. A fifth run checks C12 COV with64MiB,
the archived budget. This tests resource scaling without rerunning old code.
Archived runs used Hanhai22 EPYC9654; new runs use Hanhai25 EPYC9755. Runtime
ratios across clusters or budgets cannot isolate this kernel change.

Relevant saved controls:

| Case | Job | Steps/HVP | Application s | GNU process s | GNU RSS GiB |
|---|---|---:|---:|---:|---:|
| C12 COV |1264892|4/36|17.208014|17.40|2.5721|
| C12 ION(0:2) |1265266_1|4/44|104.286|107.11|2.824|
| C12 FULL |1265114_2|4/40|102.727|112.47|3.573|
| C14 COV |1265059|4/46|306.172876|306.62|5.7971|

Older C14 control1264476 took1190.925493s application time,4/46,energy
-535.672671742972Eh. It must not replace the better306.17s saved control.
See `benchmarks/history/c12_c14_comparison.csv` and
`results/diagnostics/c12-amplitude-20261017/README.md` for frozen evidence.

Array293257_[0-4] submits these five cases concurrently on Hanhai25.
Native results and numerical comparisons will be recorded after completion.
Build time is outside the native resource timer. Frozen input SHA256 values:

| Case | SHA256 |
|---|---|
| C12 COV |092d545498678b0cf772c5c84b9eab674b694375558cd391a1b739b7bdcc35eb|
| C12 ION(0:2) |48b44d64b773d4abaf388fcf728d478a09ab76f5426f889822b1dbc8e031432f|
| C12 FULL |3f612eab8d7d74fc6e13f1ffe1c9b96864366d664c4dadc7e96dcd4e155287b9|
| C14 COV |c9600592898c7cadf087d9cb3f77854a83eddcccc51962ac3c20bb0f393520b3|
