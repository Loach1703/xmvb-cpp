# Exact channel-support census and grouped contraction

No old molecular executable is rerun. Job293187 on Hanhai25 streams the actual
production COV support cursor for CAS(16,16), then samples every257th exchange
quotient target. It checks the prefix envelope against all supported alpha
and beta single replacements at each sampled node. CAS(8,8) is exhaustive.
This is a structural/FLOP diagnostic, not an end-to-end performance result.

| Quantity | C16 |
|---|---:|
| Spin strings per spin |12870|
| Quotient intermediate targets |50756849|
| Depth-zero quotient targets |17382365|
| Samples |197498|
| Unrestricted channels |136|
| Mean exact input channels |39.2394|
| Mean occupation-allowed output channels |105.9445|
| Mean conservative prefix envelope |60.7388|
| Input-restricted multiplication FLOPs / dense |0.288525|
| Input+output-restricted multiplication FLOPs / dense |0.224849|
| Exact input-restricted Gram FLOPs / dense |0.154729|
| Prefix-envelope multiplication FLOPs / dense |0.446609|
| Prefix-envelope Gram FLOPs / dense |0.369122|

The exact masks and both-side ratios are estimates for later algorithms, not
implemented savings. The first implementation packs rows sharing a prefix
envelope and multiplies by the corresponding columns of K. Every output
channel and the existing full Gram are retained. Generic selected/ionic
families retain the exact unrestricted operation; no covalent certificate
is applied to them. There is no floating-point screening threshold.

The existing qualified one-step control293167 (no grouped channel product)
completed in803.33s GNU-time wall,11700424KiB peak RSS,7 orbital HVPs and one
accepted step. Its energy decrease is0.87996739814116154Eh from initial
-615.12694494749928Eh. It is intentionally not converged. Do not confuse this
with intermediate293141's797.84s or phase293142's449.90s macro-only time.

New qualification and native/phase jobs will be recorded below. Kernel FLOPs
do not establish improvement: retain this change only after numeric regression
and the same C16 one-step end-to-end comparison.

Qualification293189 passes6/6 tests in6.54s. Final qualification293190 (adds
channel-work and shared-kernel-byte diagnostic counters) passes6/6 in6.49s.
Frozen executable SHA256:
`899515c153f5c4bead18fb8b3fd6841c463ba77dc75c551c39cace015a88f2cc`.
Native one-step293191 runs on anode014,32 cores,4GiB workspace, original input,
unchanged thresholds. Existing control293167 ran on anode003; both are EPYC9755.
The completed result is recorded below. The experiment remains on a work
branch; the previous qualified master is unchanged.

## Completed grouped-GEMM control

Native293191 completes its intentional one-step limit in807.78s with
11108704KiB GNU-time peak RSS. Initial/candidate CI take183.508/160.901s;
14/13 actions and7 orbital HVPs are unchanged. Accepted molecular energy is
-616.0069123456398Eh (computed from the logged initial energy and decrease).
There is **no end-to-end speedup** over803.33s control293167. The smaller
GEMM count did not justify declaring the implementation optimized.

## Direct ragged layout

The next implementation removes dense source-channel generation followed by
packing, writes scatter values directly into grouped matrices, forms smaller
channel Grams, and adds directional forcing directly into sigma when no
separate forcing output is requested. It retains full output channels and
charges the actual variable-width storage to the same4GiB action workspace.
The mathematical derivation is in `article/exterior_column_schedule.md`.

Qualification293200 passes6/6 in6.60s. Qualification293203 passes6/6 in6.82s,
additionally checking ragged planner boundaries and a C8 multi-RHS derivative
fixture large enough to exercise the per-worker grouped Gram. Frozen binary:
`26652e9d1abe9cd1adace1fefe19336dcefca17af41272513facf086be066e58`.
Native one-step293205 runs on anode014; bounded420-second profiling job293206
runs on anode031. Both use32cores. No historical version is rerun;
profiling is for attribution, not native speed claims. The one-step test does
not establish SCF convergence.

Native293205 completed in762.94s with11215844KiB GNU-time peak RSS,7 orbital
HVPs,4 keyframes and11 microiterations. The accepted energy is
-616.00691234564522Eh,4.8e-12Eh from control293167. Total time is5.03% lower
and peak RSS4.14% lower. CI/HVP/keyframe totals are343.011/189.207/55.985s;
the corresponding control totals are358.592/224.766/65.086s. The remaining
unpartitioned time is174.737s versus154.886s; it is not attributed to a
particular phase without further measurements.

Bounded capture293206 has no reported errors or warnings; its child
experiment is intentionally unclosed at the420-second timeout. Exclusive CPU
sample shares: team/barrier waiting39.87%, sigma pullback16.44%, scatter6.69%,
BLAS compute kernels6.83%. These are not wall-time fractions or a matched
full-phase before/after profile.

The next change traverses compact beta target rows and fuses exchange weights
into the channel products. It removes full-halo scans of absent target rows,
a full-image scaling pass and one team barrier. Qualification and C16 timing
will be recorded after execution; no speedup is inferred from the code change.

Qualification293253 passed6/6 tests in6.57s. At the user's request, follow-up
native/profile jobs293254/293255 were cancelled after four minutes and replaced
by C12/C14 array293257. Their incomplete timings are not performance results.
See `benchmarks/c12-c14-channels-20261018/README.md`.
