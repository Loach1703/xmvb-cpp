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
Full-step wall time and RSS are pending. Early scalar CI events agree, but
they do not establish end-to-end speedup. The experiment is kept on a work
branch pending this measurement; the previous qualified master is unchanged.
