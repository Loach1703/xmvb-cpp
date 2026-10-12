"""Summarize scalar CIAH accounting; no manual trajectory transcription."""
import csv
import json
import pathlib
import sys

run = pathlib.Path(sys.argv[1])
events = {}
for line in (run / "run.err").read_text().splitlines():
    if not line.startswith("NEO_DIAG,"):
        continue
    _, macro, _, event, payload = line.split(",", 4)
    values = dict(item.split("=", 1) for item in payload.split(",") if "=" in item)
    events.setdefault(event, []).append((int(macro), values))
norms = [v for _, v in events.get("ciah_structure_norm", [])]
ritz = [v for _, v in events.get("ciah_ritz_cost", [])]
rows = list(csv.DictReader((run / "neo.tsv").open(), delimiter="\t"))
summary = {
    "macros": len(rows),
    "hvp": sum(int(row["orbital_hvp_actions"]) for row in rows),
    "structure_norm_calls": len(norms),
    "structure_norm_seconds": sum(float(v["seconds"]) for v in norms),
    "structure_metric_actions": sum(int(v["actions"]) for v in norms),
    "structure_metric_iterations": sum(int(v["iterations"]) for v in norms),
    "max_metric_iterations": max((int(v["iterations"]) for v in norms), default=0),
    "max_certified_relative_residual": max(
        (float(v["relative_residual"]) for v in norms), default=0),
    "ritz_calls": len(ritz),
    "ritz_seconds": {key: sum(float(v[key]) for v in ritz) for key in
                     ("projection_seconds", "spectrum_image_seconds", "norm_seconds")},
    "tail": [],
}
for macro, end in events.get("ciah_micro_end", [])[-8:]:
    micro = [v for k, v in events.get("ciah_micro", []) if k == macro]
    summary["tail"].append({
        "macro": macro, "stop": end["reason"],
        "micro": len(micro), "source": float(end["source"]),
        "minimum_predicted": min(float(v["gradient"]) for v in micro),
        "last_predicted": float(end["predicted"]),
        "clipped": sum(float(v["step_scale"]) < 1 for v in micro),
    })
print(json.dumps(summary, indent=2))
