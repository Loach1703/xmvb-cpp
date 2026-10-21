#!/usr/bin/env python3
"""Bound a first-CI benchmark; stop externally after its completion event.

Run on a Slurm compute node. GNU time survives the signal to its direct
executable child and records that child's RSS. No production stopping rule
or SCF threshold is changed. A completed CI phase is not SCF convergence.
"""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import threading
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds", type=float, required=True)
    parser.add_argument("--workspace-mib", type=int, default=4096)
    parser.add_argument("binary")
    parser.add_argument("input")
    args = parser.parse_args()
    result = {"scope": "first_ci", "state": "running", "workspace_mib": args.workspace_mib}
    lock = threading.Lock()
    start = time.monotonic()
    with open("output.out", "wb") as output, open("run.err", "wb") as errors:
        process = subprocess.Popen(
            ["/usr/bin/time", "-v", "-o", "resource.time", "stdbuf", "-oL", "-eL",
             args.binary, args.input, "--optimizer-backend", "ciah",
             "--ci-workspace-mib", str(args.workspace_mib)],
            stdout=output, stderr=subprocess.PIPE, start_new_session=True)

        def stop(reason):
            with lock:
                if result["state"] != "running":
                    return
                result.update(state=reason, wall_seconds=time.monotonic() - start)
            # time's one direct child is stdbuf, which execs the executable.
            # Keep time alive to collect the native child's wait4 usage.
            children = Path(f"/proc/{process.pid}/task/{process.pid}/children")
            try:
                for child in children.read_text().split():
                    os.kill(int(child), signal.SIGTERM)
            except (FileNotFoundError, ProcessLookupError):
                pass

        timer = threading.Timer(args.seconds, stop, args=("diagnostic_timeout",))
        timer.start()
        try:
            for line in process.stderr:
                errors.write(line)
                errors.flush()
                elapsed = time.monotonic() - start
                if b",ciah_ci_iteration," in line:
                    result.setdefault("first_ritz_wall_seconds", elapsed)
                    result["last_ritz"] = line.decode().strip()
                if b",ciah_ci_solved," in line:
                    result["ci_completion"] = line.decode().strip()
                    stop("initial_ci_completed")
            code = process.wait(timeout=10)
        finally:
            timer.cancel()
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        if result["state"] == "running":
            result.update(state="program_exited", wall_seconds=time.monotonic() - start)
        result["wrapper_exit"] = code
        result["capture_wall_seconds"] = time.monotonic() - start
        Path("result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)
    return 0 if result["state"] in {"initial_ci_completed", "diagnostic_timeout"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
