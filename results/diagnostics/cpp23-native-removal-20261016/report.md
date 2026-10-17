# Remove native L-BFGS and require C++23

Parent: `5f08336b`. Remove ISCF=5 rather than redirecting it to another
algorithm. P-LBFGS remains the default, ISCF=8; TNHVP=7 and CIAH=9 are unchanged.

## Cleanup

- Delete the raw-coordinate inverse recursion, its 100-pair history and
  original L-BFGS backend. Rename the retained backend to `backends/plbfgs.*`.
- Delete the exclusive Moré--Thuente search, interpolation helper and result
  type. Keep the shared Armijo search used by the supported optimizers.
- Remove the input/driver enum entries, CLI selection, dispatcher branch,
  native positive test and raw-only iteration-report branches.
- ISCF=5 and CLI `lbfgs` now fail explicitly. Regression tests verify rejection
  rather than accepting a compatibility alias. The three versioned decks
  previously specifying `ISCF=5` now explicitly request `ISCF=8`; original result archives
  are not modified.

Source comparison confirms that the complete retained `run_plbfgs_backend`
implementation and retained Armijo/retraction implementation are byte-identical
to the parent. No physical derivative, preconditioner, strategy or convergence
tolerance changed. Generic L-BFGS/secant operations still needed by P-LBFGS
and TNHVP are not the deleted native backend and remain available.

## Toolchain

Require C++23 and CMake 3.20 or newer, retaining `CXX_STANDARD_REQUIRED=ON`,
extensions OFF and no `-fpermissive`. Update both production and standalone
test compile requirements. GCC 14.1 on Hanhai22 was used for the full build.
The saved numerical-library flags contain `-std=c++23`, `-O3 -DNDEBUG`, `-g1`
and `-fno-omit-frame-pointer`; no numerical flags were relaxed.

C++20 was the preceding compiler-coverage choice, not an algorithmic
requirement. The current toolchain supports enabling C++23. Specific newer
facilities should be chosen according to their actual implementation status,
not assumed universally present just because the language mode is enabled.
See [GCC14](https://gcc.gnu.org/gcc-14/changes.html) and
[CMake3.20 CXX_STANDARD](https://cmake.org/cmake/help/v3.20/prop_tgt/CXX_STANDARD.html).

## Verification

Hanhai22 job **1264229**, anode23: full build and **34/34 regressions passed**.
This includes actual P-LBFGS F2 convergence, CIAH HAO/OEO and matched-seed
reference, state averaging, keyword/CLI rejection, trace/text/orbital output,
common solvers, CASCI, transport and linear algebra. Both developer diagnostics
compile. No production performance claim is made from these correctness tests.

Source hashes match 497 local files. The two differences are unchanged,
pre-existing local RI edits excluded from this commit and cluster transfer.
Raw build/test logs, flags, dependency cache, source/binary hashes and library
linkage symbols are retained in `evidence.tar.gz`; accounting is separate.

Executable SHA-256:
`ded4a97f304e98620351bb8ebd7a974d3c61517eaf4295ee765fac59d5d7a599`.
