# Same-guess orthogonal CASSCF validation

`export_fixture.py` creates one all-electron Cartesian AO integral and complete
MO-frame fixture from a documented orbital seed. `pyscf_newton.py` and the standalone
`casscf_fixture` executable read that fixture. No solver regenerates
AO integrals or a new orbital guess. The primary reference is PySCF's coupled
second-order `newton_casscf` (`mcscf.CASSCF(...).newton()`), not the default
one-step orbital-only AH method. PySCF documents `gen_g_hop` as a joint
orbital/CI action and its kernel records `tot_hop` in the macro callback.
For a separate optimized-standard comparison, `pyscf_newton.py --solver ah`
retains the ordinary `mcscf.CASSCF` orbital-only AH update. The default
`--solver newton` keeps the coupled reference. These methods must be labeled
separately; their HVP actions do not have the same meaning. Both read the
identical fixture and report solver wall time, callback trajectory, final
energy, and final spin expectation. Missing AH `tot_hop` remains JSON null.

Cases are grounded in `testdata/vbscf/F2_OEO_PYSCF_LOW.xmi` and
`testdata/vbscf/FeCl2.xmi`. F2 is neutral, 1.4 Å, cc-pVDZ, CAS(2e,2o),
singlet, with 8 inactive orbitals. FeCl2 is neutral and linear at Fe–Cl
2.136319 Å, def2-TZVP, CAS(10e,8o), quintet, with 25 inactive orbitals.
For F2, the complete RHF frame is reordered to put the highest occupied `Ag`
sigma bonding orbital and lowest virtual `B1u` sigma antibonding orbital in
the active window. A separate D2h molecule labels these two columns, while
the fixture molecule and CASSCF calculations keep symmetry disabled. The
selected irrep weights and original MO indices are recorded. This avoids an
accidental frontier-pi CAS(2,2). For FeCl2, the first 33 orbitals come from
the archived full-AO `$GUS` input; the qualified XMVB `input.molden` preserves
the AO shell order. The `fecl2_seed.py` QR procedure projects active columns
off the 25 core columns and completes the virtual complement. FeCl2 does not
run a preparatory SCF calculation. Its recorded single-determinant energy is
only a property of that seed, not a converged HF result.
The archived Molden atom labels are indexed (`Fe1`, `Cl2`, `Cl3`); the exporter
validates their nuclear charges, Cartesian 136-AO layout, and exact geometry.
The bounded layout check `inspect_fecl2_layout.sbatch` measured a 2.28e-10
maximum overlap difference against the source GBS basis on Hanhai25 job 291879.

`ao-mo.bin` is little-endian `XMCASAO1`: 8-byte magic; five uint32 values
(`nbf`, `ncore`, `ncas`, active electrons, `2S`); nuclear repulsion as float64;
then row-major float64 `S`, `hcore`, and full `U`, each `nbf x nbf`; finally
float64 `int2e` in PySCF s8 lower pair-of-pairs order. Pair `(i,j)` uses
`i*(i+1)/2+j` for `i>=j`, matching `core::PairIndex`. The ERI payload has
`n_pair*(n_pair+1)/2` values. There is no dense four-index AO array in the
fixture or the C++ reader. The reader builds the shared `AoPairGraph` once.
The JSON sidecar retains source hashes, Molden layout hash for FeCl2, PySCF
version, MO Gram error, and the independent CASCI energy at the exact saved
initial frame.

The `--initial-only` CLI and `check_initial.py` enforce initial CASCI parity
within 1e-8 Eh before either full optimization. A different
initial energy means integral, basis, spin, or orbital-frame mismatch and makes
the optimization comparison invalid. After both solvers converge, compare
energies within 1e-6 Eh. Macro/HVP counts and wall/RSS are recorded separately:
the two solvers have different update policies, so equal counts are not
expected. `run.sbatch` records executable/input hashes, modules, compiler,
affinity, GNU time resource reports, native trajectories, and raw output.
It takes `SOURCE_DIR CASE` (for example, the qualified archive mirror and
`F2`), verifies `parent-revision.txt` and `source.sha256`, and uses the shared
`SOURCE_DIR/build/src/casscf_fixture` built by the separate qualification gate.
It does not configure or build a second copy per case. The PySCF runtime uses
Anaconda 25.06 Python 3.13 with the archived PySCF 2.14 wheel overlay from
the existing Hanhai25 benchmark; the wheel hash is retained. The submitter
should also capture `sacct` for the application steps and verify the live
partition/QOS/limits before submission. Builds and molecular workloads must
execute on a Hanhai25 Slurm compute node.

After a shared `fixture/` already exists, `reference.sbatch SOURCE_DIR
FIXTURE_DIR CASE` runs only the standard AH PySCF solver, with the same 32-core
binding and wheel. It does not repeat integral export, preparatory SCF, native
CASSCF, or the coupled Newton reference. The expected case is checked against
the saved fixture; output lives in a new `casscf-ah-CASE-JOBID` directory.
