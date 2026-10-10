"""Recompile the two edited kernels for GCC vectorization evidence, not timings."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

if not os.environ.get("SLURM_JOB_ID"):
    raise SystemExit("Run on a Slurm compute node")
commands = json.loads(Path(sys.argv[1]).read_text())
root = Path(sys.argv[2]).resolve()
for name in ("sigma.cpp", "direct_shell.cpp"):
    matches = [entry for entry in commands if entry["file"].endswith(
        "/orthogonal_ci/sigma.cpp" if name == "sigma.cpp" else "/libcint/direct_shell.cpp")]
    if len(matches) != 1:
        raise SystemExit(f"Expected one compile command for {name}")
    entry = matches[0]
    command = shlex.split(entry["command"])
    command[command.index("-o") + 1] = str(root / f"{name}.o")
    # This diagnostic object is not linked into the measured executable.
    command.append(f"-fopt-info-vec-optimized={root / (name + '.vector.txt')}")
    with (root / f"{name}.compiler.txt").open("w") as log:
        subprocess.run(command, cwd=entry["directory"], stdout=log,
                       stderr=subprocess.STDOUT, check=True)
