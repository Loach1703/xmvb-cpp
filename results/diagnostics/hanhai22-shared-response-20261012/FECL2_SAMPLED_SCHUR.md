# FeCl2: sampled Schur softening and missing coarse samples

## Scope

Hanhai22 1249478, anode31, 32 bound cores, original FeCl2 orbitals,
exact integrals and analytic diagonal. Numerical policy is unchanged.
All five unit executables passed. Native convergence: 14 macros / 281 HVPs,
-2181.624464869630 Eh, GNU wall 9.93 s, peak process RSS 492.59 MiB.
The endpoint is 1.4103e-8 Eh above the PySCF coupled-Newton reference.
No old benchmark was rerun, and no iron long run was launched.

## Diagnostic model

The existing coarse cache already owns independent orbital and structure
samples and both their complete Hessian and physical metric images. Split
these into columns V_o and V_s, then form small sampled blocks

$$
A_V=V_o^THV_o,\quad C_V=V_s^THV_s,\quad
B_V=V_s^THV_o,\quad M_o=V_o^TGV_o,\quad M_s=V_s^TGV_s.
$$

Here H and G are abstract coupled operators, not assembled production
matrices. When sampled C_V is positive, eliminate only its sampled space:

$$
S_V=A_V-B_V^TC_V^{-1}B_V.
$$

All spectra use their respective physical metric. The C_V solve is
residual-certified; no eigenvalue clipping or regularization is introduced.
The curvature removed is measured by the largest generalized eigenvalue of
(B_V^TC_V^{-1}B_V,A_V) when A_V is positive. None of these quantities changes
a step, tolerance or stopping condition. Storage is at most the existing
60-column coarse model, not a complete orbital or structure Hessian.
No additional operator or integral evaluation occurs.

Independent unit tests use an explicit four-dimensional coupled problem,
nonidentity metric, analytically known Schur complement and known removal
fraction. They also check that negative structure curvature is reported,
not shifted away. Sampled positivity does not certify full-space positivity;
sampled spectra are not full-space condition numbers.

## Results

| Macro | Orbital / structure samples | Structure spectrum | Orbital spectrum | Sampled relaxed spectrum | Maximum fraction removed |
|---|---:|---:|---:|---:|---:|
| 10 | 21 / 21 | 1.038–5.330 | 0.01859–314.5 | 0.009113–143.9 | 99.252% |
| 11 | 22 / 22 | 1.019–4.968 | 0.01739–605.7 | 0.007283–270.0 | 99.574% |
| 12 | 23 / 23 | 1.002–5.363 | 0.01588–845.4 | 0.008312–306.3 | 99.232% |
| 13 | 27 / 27 | 0.980–5.058 | 0.01705–996.4 | 0.01036–262.6 | 98.005% |
| 14 | 29 / 18 | 0.992–4.849 | 0.01666–1042.1 | 0.01112–384.3 | 98.824% |

The structure samples have a moderate physical spectral spread, whereas
coupling nearly cancels orbital curvature in some sampled directions.
Macros 1–4 have negative sampled relaxed curvature; the diagnostic retains
it. An independent diagonal model cannot represent this cancellation.
The existing coupled coarse inverse does represent it inside its admitted
sample space; this observation alone does not justify another Woodbury
coupling approximation (the earlier such trial regressed).

## Concrete sample-admission issue

The rank certificate rejects nine samples in macro 8 and eleven in macro
14. Their remaining Euclidean direction fractions are respectively
0.0873–0.6621 and 0.2055–0.4809. Thus these are not locally negligible
innovations. Rejection follows the *global original-source Gram* certificate:
the operator images were evaluated before recursively transforming nearly
dependent sources. Large local innovation does not certify that transformed
cached images are accurate, so the existing rejection must not simply be
weakened or removed.

The next mathematically distinct repair to test is component-wise operator
reconstruction with freshly evaluated innovations. For each orbital or
structure component x and an already orthonormal cached basis V,

$$
a=V^Tx,\qquad z=x-Va,\qquad Hx=(HV)a+Hz.
$$

Evaluate the operator on the independent normalized innovation *before*
admitting its image, instead of obtaining its image by subtracting old
already-evaluated near-dependent source images. Apply the identical
decomposition to metric images. This can replace, not supplement, the
existing orbital/structure operator evaluation. Linearity preserves the
complete AH action; no coupling is discarded and no rank threshold changes.
Dependent residual components must still be accounted for, not truncated.

Such a fresh sample avoids amplification by the old original-source QR
factor because its operator image was evaluated directly. Exact cached
operator reconstruction, symmetry and the declared operator accuracy must
be verified independently before any speedup claim. This algorithm is not
implemented by this diagnostic commit; the next action is a bounded test
of reconstruction/admission, not an increase of the AH budget.

Raw logs, hashes and GNU resource records: fecl2-ritz-cost-1249478/native/.
summarize_ritz_cost.py extracts the sampled spectra programmatically.

## Fresh-innovation reconstruction: tested and rejected

The proposed component reconstruction was implemented. Orbital and structure
components were separately projected against cached orthonormal samples,
their independent normalized residuals were evaluated directly, and full
operator images were reconstructed by linearity. Small unadmitted residuals
were still evaluated; no residual component was truncated. Only directly
evaluated independent images were admitted. Complete metric images used
the same fresh directions. Original rank criteria, action precision,
CIAH policies and nonlinear acceptance were unchanged.

Independent dense unit tests checked cross-block coupling, the chained-source
case rejected by the original Gram certificate, full operator reconstruction,
one action per nonzero innovation and evaluation of a below-admission
residual with a large Hessian image. All five unit executables passed in
both experiment jobs.

| FeCl2 path | Job | Macro / orbital HVP | Energy / Eh | GNU wall / s | RSS / MiB | Expansion scope / s |
|---|---:|---:|---:|---:|---:|---:|
| Archived retained diagnostic | 1249478 | 14 / 281 | -2181.624464869630 | 9.93 | 492.6 | 2.447 |
| Fresh orbital and structure components | 1249491 | 14 / 280 | -2181.624464867051 | 10.19 | 511.6 | 2.837 |
| Fresh structure component only | 1249494 | 14 / 281 | -2181.624464866901 | 10.05 | 500.0 | 2.755 |

All observations are native runs on Hanhai22 anode31 at 32 bound cores;
the old binary was not rerun. All endpoints agree with the PySCF
coupled-Newton reference to about 1.7e-8 Eh. The all-fresh path ends with
28 orbital and 28 structure samples, resolving the old 29/18 admission
imbalance. It nevertheless does not reduce macro steps or materially
reduce HVP counts. Additional metric evaluations and reconstruction work
coexist with a larger measured expansion scope; that scope alone does not
attribute its increase to a specific instruction or allocation.

This rejects sample admission as the dominant FeCl2 tail bottleneck.
The prior rejection is a legitimate cached-image conditioning certificate,
not evidence that the physical HVP loses those directions: rejected samples
were absent from the *preconditioner*, not from the complete coupled AH
action. The new path made richer samples trustworthy but did not demonstrate
an effective convergence improvement. Its method, driver refactor and trial
unit test were deleted locally and from the cluster source mirror. There is
no retained switch or fallback. The mirror executable was restored from
frozen 1249478 without running another old benchmark.

The next unresolved question is representation consistency of the approximate
inverse: the analytic orbital rotation diagonal is built in an accepted
orthogonal carrier frame, whereas the complete coupled Hessian acts on raw
nonorthogonal orbital/structure coordinates. Existing source/frame transforms
are present and independently tested, so their absence must not be asserted.
What remains to inspect is whether the carrier connection belongs in a
coupled congruence of the *preconditioning model*. This is not a proposed
change to orbital retraction, the exact Hessian or stopping rules. The
connection audit already exists, but a production congruence has not been
implemented or validated. No new iron long run was launched.
