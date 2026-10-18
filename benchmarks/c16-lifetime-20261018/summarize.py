#!/usr/bin/env python3
"""Parse saved first-CI JSON, GNU time and first complete sigma-sweep records."""
import csv
import hashlib
import json
from pathlib import Path
import re
import sys


def fields(line):
    return dict(re.findall(r"(\w+)=([^,\s]+)", line))


def summarize(run):
    budget_run=run.name.startswith("mib-")
    root=run.parent if budget_run else run
    native=run if budget_run else run/"native"
    phase_path=root/("phases-4096" if budget_run else "phases")/"run.err"
    result = json.loads((native / "result.json").read_text())
    completed = fields(result.get("ci_completion", ""))
    usage = (native / "resource.time").read_text()
    peak = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", usage)
    lines = phase_path.read_text().splitlines() if phase_path.exists() and (not budget_run or run.name=="mib-4096") else []
    sigma = next((fields(x) for x in lines if x.startswith("XMVB_TILE,")), {})
    columns = next((fields(x) for x in lines if x.startswith("XMVB_PRODUCT_ACTION,")), {})
    release = next((fields(x) for x in lines if x.startswith("XMVB_COLUMN_RELEASE,")), {})
    slots = int(release.get("factor_slots", "0"))
    return {
        "job": root.name.removeprefix("run-").removeprefix("budgets-"),
        "workspace_mib": result.get("workspace_mib",256),
        "node": (root / "node.txt").read_text().strip(),
        "threads": 32,
        "scope": result["scope"],
        "state": result["state"],
        "endpoint_wall_seconds": result["wall_seconds"],
        "first_ritz_wall_seconds": result.get("first_ritz_wall_seconds", ""),
        "ci_solver_seconds": completed.get("seconds", ""),
        "ci_iterations": completed.get("iterations", ""),
        "ci_block_actions": completed.get("actions", ""),
        "active_ci_energy": completed.get("energy", ""),
        "peak_rss_kib": int(peak.group(1)) if peak else "",
        "first_sigma_stage_sum_seconds": sum(float(sigma[k]) for k in ("prepare", "read", "compute", "write")) if sigma else "",
        "first_sigma_targets": sigma.get("targets", ""),
        "input_evictions": columns.get("input_evictions", ""),
        "output_evictions": columns.get("output_evictions", ""),
        "factor_density": int(release["factor_entries"]) / slots if slots else "",
        "wrapper_exit": result["wrapper_exit"],
        "binary_sha256": (root / ("binary.sha256" if budget_run else "hashes.txt")).read_text().split()[0],
        "manifest_sha256": hashlib.sha256((root / "source.sha256").read_bytes()).hexdigest(),
    }


rows = [summarize(Path(path)) for path in sys.argv[1:]]
if not rows:
    raise SystemExit("usage: summarize.py RUN_DIR [RUN_DIR ...]")
writer = csv.DictWriter(sys.stdout, fieldnames=rows[0])
writer.writeheader()
writer.writerows(rows)
