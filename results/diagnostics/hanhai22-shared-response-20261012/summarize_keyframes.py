#!/usr/bin/env python3
"""Read saved scalar traces; no optimization jobs or counterfactual claims."""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

run = Path(sys.argv[1])
events = defaultdict(list)
for line in (run / "run.err").read_text().splitlines():
    if line.startswith("NEO_DIAG,"):
        _, macro, attempt, event, payload = line.split(",", 4)
        values = dict(item.split("=", 1) for item in payload.split(",") if "=" in item)
        events[event].append(dict(macro=int(macro), attempt=int(attempt), **values))
rows = list(csv.DictReader((run / "neo.tsv").open(), delimiter="\t"))
costs = events["ciah_keyframe_cost"]
active = events["ciah_keyframe_active_cost"]
summary = {
    "run": str(run),
    "macros": len(rows),
    "hvp": sum(int(r["orbital_hvp_actions"]) for r in rows),
    "keyframes": sum(int(r["keyframes"]) for r in rows),
    "keyframe_calls_timed": len(costs),
    "keyframe_seconds": sum(float(v["seconds"]) for v in costs) if costs else None,
    # Nested timers: active is part of keyframe, not an additional wall cost.
    "keyframe_active_seconds": {
        key: sum(float(v[key]) for v in active) if active else None
        for key in ("seconds", "setup", "action", "adjoint")},
    "keyframe_hs_actions": sum(int(v["actions"]) for v in active) if active else None,
    "keyframe_max_diagonal_bytes": max(
        (int(v["diagonal_bytes"]) for v in active), default=None),
    "keyframe_decisions": events["ciah_keyframe"],
    "micro_ends": events["ciah_micro_end"],
    "macro_trials": events["ciah_macro_trial"],
    "limits": "Observed corrections do not establish the outcome of skipping a keyframe. "
              "Scalar-audit wall time is not the native benchmark time.",
}
print(json.dumps(summary, indent=2))
