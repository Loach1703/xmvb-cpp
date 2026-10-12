# TNVB active-space adjoints and independent backend

This is the mathematical boundary for the fixed-bond SU(2), DIRECT, FULL
backend. This interface was derived before implementation; numerical evidence
and acceptance evidence are tracked in `tnvb_acceptance_report.md`. No determinant
response object belongs to this backend.

## Energy and structure coordinates

The C++ orbital/integral preparation supplies a positive definite active
overlap S, effective one-electron integrals h (including the inactive field),
chemist's integrals g_ijkl=(ij|kl), and the inactive reference E11. The worker
Hamiltonian has zero core energy. Its generalized eigenvalue is Eact, so

    Eelectronic = E11 + Eact
    Etotal = Enuclear + E11 + Eact.

The worker constructs individually normalized structure MPSs and solves
H c = Eact Sstructure c, with c^T Sstructure c = 1. H and Sstructure are
reported in these same coordinates. The DIRECT total MPS uses exactly c.
No additional structure normalization is applied during the C++ pullback.

The C++ paired-CSF convention is invariant under exchanging the members of
one singlet pair. The worker adapter sorts the two indices within each active
pair before calling the reused oriented-pair core. Pair order and inactive
occupations are retained; the caller's input is not mutated. This ensures the
reported H/S rows and coefficients use the same phase convention as the C++
determinant expansion. No determinant expansion is needed for this mapping.

For an exact stationary root, changing a structure's normalization multiplies
both its row/column coordinates. The resulting derivative is proportional to
(H-Eact Sstructure)c and vanishes. Thus the envelope derivative can be computed
with the instantaneous normalized-basis coefficients held fixed; it must not
include a second normalization correction.

## Spin-free densities and permutations

In the worker's orthonormal basis, with S=R^T R and Q=R^-1,

    b_i^dagger = sum_a R_ai a_a^dagger
    horth = Q^T h Q
    gorth_abcd = sum_ijkl Q_ia Q_jb Q_kc Q_ld g_ijkl.

Block2's real spin-summed densities are

    gamma_ab = sum_sigma <a_a,sigma^dagger a_b,sigma>
    Gamma_abcd = sum_sigma,tau
                 <a_a,sigma^dagger a_b,tau^dagger a_c,tau a_d,sigma>.

Each density is divided by the total MPS norm. In particular,

    Eact = sum_ab horth_ab gamma_ab
         + 1/2 sum_abcd gorth_abcd Gamma_acdb.

The worker transforms every free index contravariantly with Q and exports

    P_ij       = gamma_nonorth_ji
    Pi_ijkl    = Gamma_nonorth_lkij
    C_dc       = sum_abef gorth_faeb Gamma3_abcdef
    D3G        = Q C Q^T.

These are column-major boundary buffers. Pi's flat index is
`i + K*(j + K*(k + K*l))`. FULL Gamma3 is contracted and released inside the
worker; it is never sent to the C++ process.

## Three active-space adjoints

With the Frobenius convention dE = sum_ij A_ij dX_ij, the real, symmetric
matrix adjoints are

    Ah_ij = P_ij
    As_ij = sym[ sum_uv h_uv Pi_v i u j + 1/2 D3G_ij - Eact P_ij ]
    Ag[p(ij),p(kl)] = 1/2 sum_(a,b,c,d in the packing orbit)
                                     Pi_d b c a.

`sym(M)=(M+M^T)/2`. As is the derivative at fixed h and g, not the
derivative of horth at a fixed MPS. Independent symmetric off-diagonal S/h
perturbations must use both matrix entries in their directional dot product.
The packed Ag is a sum over every ordered quadruple mapping to the same
eightfold-symmetric integral; repeated-index permutations appear only once.
There is no extra universal factor of 2, 4, or 8 after this accumulation.

### Overlap-adjoint derivation

Let X=dR Q. Changing the orbitals in a fixed structure polynomial gives
d|Psi> = sum_ab X_ab E_ab |Psi>, where E_ab is the spin-free excitation.
Consequently its normalized state contribution is

    2 sum_ab X_ab [ <E_ba H> - Eact gamma_ab ].

Normal ordering gives

    <E_ab H1> = sum_q horth_bq gamma_aq
               + sum_pq horth_pq Gamma_apqb
    <E_ab H2> = 1/2 C_ba + sum_uvx gorth_uvbx Gamma_auv x.

The simultaneous change of horth and gorth at fixed nonorthogonal h,g cancels
the first term of the H1 transition and the lower-rank term of the H2
transition. The remaining derivative is

    dEact = sum_ab X_ab [2 sum_pq horth_pq Gamma_b p q a
                        + C_ab - 2 Eact gamma_ab].

Since dS=R^T(X+X^T)R, transforming the free indices with Q yields As above.
The -Eact P term is the nonorthogonal norm derivative; E11 and nuclear
repulsion do not belong in it.

At finite MPS bond dimensions these formulas differentiate the represented
state under the exact-polynomial convention. They do not differentiate every
truncation and fit operation. Eact remains the H/S eigenvalue used by the
optimizer and output; the fitted state's expectation discrepancy is reported
as fit_dE. A value above 1e-7 Eh is a production warning, not an optimizer stop
or convergence condition. Formal numerical verification separately requires
fit_dE<=1e-7; this distinction does not change production behavior.

## Shared C++ pullback and RI

Pass As, Ah and the packed Ag to the existing C++ active-integral pullback.
It adds the E11 inactive-density/reference terms exactly once and propagates
through orbital normalization, inactive projection, and the HAO/OEO chart.
An explicit selected structure space must leave complete_active_space=false.

For RI, g_PQ = sum_A L_AP L_AQ is reconstructed only in active space for
the worker. The C++ result retains its RI representation. Reverse propagation
uses dL = L (Ag_pair + Ag_pair^T), with the packed-pair adjoint expanded using
the existing C++ convention, and then uses the existing AO RI pullback. No
dense AO four-index tensor is reconstructed.

## Required numerical gates

Validate integral directions in S, h and packed g separately before orbital
optimization. Validate active and inactive orbital directions over several
finite-difference step sizes. Compare normalized H/S and phase-aligned
coefficients/weights with the determinant backend. HAO/OEO and exact/RI are
separate required configurations. Nonfinite data, invalid structure metrics,
worker failures, or RDM identity errors above 1e-8 remain fatal.
