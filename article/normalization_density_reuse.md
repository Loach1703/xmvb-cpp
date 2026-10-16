# Reusing the integral cotangent for normalization terms

Status: the density identity and its derivative were independently validated.
The complete production rewrite was tested in Hanhai22 jobs 1264640 and
1264642, but not retained: fewer generator contractions did not yield a faster
complete C12 run. The current production cache still stores Y and (H-E)Y.
The formulas below remain a validated option for a future larger fusion.

For real spin carriers define the oriented generator moment

$$
D(L,R)_{pq}=\langle L|E_{pq}|R\rangle,\qquad
E_{pq}=\sum_\sigma a^\dagger_{p\sigma}a_{q\sigma}.
$$

The symmetric full-matrix one-electron integral cotangent is

$$
J(L,R)=\tfrac12\left[D(L,R)+D(R,L)\right].
$$

This is already returned by the direct-CI integral adjoint. The Hamiltonian
channel embedding uses a fixed-electron sector, where the total number operator
equals the electron count times the identity, so its one-electron pullback is
exactly this transition density. In particular, $J(Y,Y)=D(Y,Y)$.

For one weighted state and a fixed Lagrange multiplier $E$, the representation
generator cotangent is

$$
G=2w\operatorname{up}\left[D(HY,Y)-EJ(Y,Y)\right].
$$

For $Z=\dot Y$ and a simultaneous physical integral direction $\dot H$,

$$
\dot G=2w\operatorname{up}\left[
D(HZ+\dot H Y,Y)+D(HY,Z)-2EJ(Z,Y)\right].
$$

The existing directional integral adjoint returns $J(Z,Y)$ in the same sigma
sweep. Its weighted one-electron response is $2wJ(Z,Y)$. Therefore the energy
normalization term in $\dot G$ is obtained by multiplying that existing small
matrix by $-2E$, without another generator contraction or carrier traversal.
The accepted normalization term similarly reuses the already computed
$wJ(Y,Y)$. This identity holds for arbitrary, nonstationary carriers and any
fixed E; it is not a residual-based switch or an approximation at stationarity.

## Accepted image representation

Retain Y and HY, rather than Y and (H-E)Y, in the existing two-component compact
cache. Both use exactly the same support and storage budget. The accepted HY
supplies the second generator moment above directly, while the cached integral
cotangent supplies all normalization terms. Left directional projection also
uses HY directly, without unpacking Y or reconstructing HY from a residual.

Nonresident primal/tangent sweeps already supply both Hamiltonian images and
therefore need no final normalization-generator traversal. If neither images
nor diagnostics are requested, their now-unused final carrier read is omitted.
Selected-space projection, physical CI response, off-diagonal state multipliers,
finite transport, and the physical gradient/HVP pullbacks are unchanged.
No closure of the selected VB subspace is required. Equal and unequal spin
sectors use the same spin-summed density identity; spin parity is not its premise.
