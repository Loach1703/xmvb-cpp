#!/usr/bin/env python3
"""ENERGY, overrides, capabilities, fit warnings and analysis/export tests."""
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def main():
    executable, fixture = sys.argv[1:]
    original = Path(fixture).read_text()
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", GOTO_NUM_THREADS="1",
               XMVB_BLOCK2_WORKERS="1", SLURM_CPUS_PER_TASK="2",
               XMVB_BLOCK2_SUM_ENERGY_TOL="1e-99", XMVB_BLOCK2_RDM_CHECK_TOL="1e2")
    with tempfile.TemporaryDirectory(prefix="cpp-tnvb-cli-") as directory:
        path = Path(directory) / "input.xmi"

        def run(configuration, iterations, options=(), forbidden=True, molden=False, iscf=5):
            text = original.replace("tnvb\ntnbond=128\ntnsumbond=256\n", configuration + "\n")
            text = text.replace("itmax=0", f"itmax={iterations}")
            text = text.replace("iscf=5", f"iscf={iscf}")
            if not molden:
                text = re.sub(r"(?mi)^molden\s*$", "", text)
            artifact = path.with_suffix(".molden")
            if artifact.exists():
                artifact.unlink()
            path.write_text(text)
            child_env = env.copy()
            if forbidden:
                child_env["XMVB_BLOCK2_FORBID_3PDM"] = "1"
            else:
                child_env.pop("XMVB_BLOCK2_FORBID_3PDM", None)
            result = subprocess.run([executable, str(path), *options], env=child_env,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if molden and result.returncode in (0, 2):
                assert artifact.is_file(), result.stdout
                exported = artifact.read_text()
                for section in ("[Molden Format]", "[Atoms]", "[GTO]", "[MO]"):
                    assert section in exported, section
                # This fixture has 18 inactive and 6 active orbitals. The existing
                # export must retain all of them, including in ENERGY mode.
                assert len(re.findall(r"(?m)^Sym=", exported)) == 24, exported
                assert "Molden output" in result.stdout, result.stdout
            return result.returncode, result.stdout

        for configuration, itmax, options in [
                ("TNVB TNBOND=128", 0, ()),
                ("TNVB TNBOND=128", 200, ("--max-iterations", "0")),
                ("TNVB TNBOND=128 EIGENSOLVER=DAVIDSON", 0, ("--eigensolver", "dense"))]:
            status, output = run(configuration, itmax, options, molden=True)
            assert status == 0 and "Termination reason                 : energy_only" in output, output
            for label in ["OVERLAP OF VB STRUCTURES", "HAMILTONIAN OF VB STRUCTURES",
                          "COEFFICIENTS OF STRUCTURES", "WEIGHTS OF STRUCTURES",
                          "ORBITALS IN PRIMITIVE BASIS FUNCTIONS", "Fixed-orbital ENERGY"]:
                assert label in output, label
            assert "[TNVB-RDM]" not in output and "[TNVB-SUM]" not in output, output
            assert "VBSCF converged" not in output and "Final gradient |g|" not in output, output
            assert "COMPUTED NATURAL ORBITALS" not in output, output
            assert "XMVB ATOMIC POPULATION ANALYSIS" not in output, output
        for configuration, options, message in [
                ("TNVB TNBOND=128", ("--max-iterations", "1"), "requires explicit TNSUMBOND"),
                ("TNVB TNBOND=128 TNSUMBOND=64", ("--max-iterations", "1"), "requires explicit TNSUMBOND"),
                ("TNVB TNBOND=128", ("--eigensolver", "davidson"), "dense H/S"),
                ("TNVB TNBOND=128", ("--optimizer-backend", "neo"), "ISCF=5"),
                ("TNVB TNBOND=128", ("--optimizer-backend", "nonredundant_truncated_newton"), "ISCF=5"),
                ("TNVB TNBOND=128", ("--max-iterations", "-1"), "ISCF=5")]:
            status, output = run(configuration, 0, options)
            assert status != 0 and "XMVB-CPP error:" in output and message in output, output
            assert "[TNVB-MPS]" not in output, "capability failure launched tensor calculations"
        # Each structure is exact here, but cap 16 truncates the sum (exact rank 21).
        status, output = run("TNVB TNBOND=16 TNSUMBOND=16", 1, forbidden=False)
        assert status in (0, 2) and "XMVB-CPP error:" not in output, output
        assert "[TNVB-WARNING] fit_dE=" in output and "continuing RDM/gradient" in output, output
        assert "[TNVB-RDM]" in output and "Total Energy:" in output, output
        # A completed full-gradient evaluation retains existing postprocessing
        # even when the accepted-step cap stops optimization before convergence.
        status, output = run("TNVB TNBOND=128 TNSUMBOND=256", 1,
                             forbidden=False, molden=True)
        assert status in (0, 2) and "XMVB-CPP error:" not in output, output
        for label in ("COMPUTED NATURAL ORBITALS", "XMVB ATOMIC POPULATION ANALYSIS",
                      "POPULATION AND CHARGE", "MULL.POP.", "LOW.POP."):
            assert label in output, label
        # Input ISCF=8 and its CLI override use the existing block optimizer.
        # ENERGY still must not build a total MPS or call any RDM.
        status, output = run("TNVB TNBOND=128", 0, iscf=8)
        assert status == 0 and "energy_only" in output, output
        assert "[TNVB-RDM]" not in output and "[TNVB-SUM]" not in output, output
        for iscf, options in [(8, ()), (5, ("--optimizer-backend", "block_lbfgs"))]:
            status, output = run("TNVB TNBOND=128 TNSUMBOND=256", 1, options,
                                 forbidden=False, molden=True, iscf=iscf)
            assert status in (0, 2) and "XMVB-CPP error:" not in output, output
            for label in ("block-LBFGS (nonredundant orbital-block inverse)",
                          "projected |g|_inf", "Final projected |g|_inf",
                          "COMPUTED NATURAL ORBITALS", "XMVB ATOMIC POPULATION ANALYSIS"):
                assert label in output, label
            assert "block_lbfgs_line_search_failed" not in output, output
    print("TNVB CLI ISCF=5/8, ENERGY, override, capability, fit-warning and analysis/export tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
