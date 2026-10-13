"""Summarize existing scalar diagnostics; does not run chemistry workloads.

Rows are named timer scopes, not an additive wall-time decomposition.
Usage: python3 summarize-remaining-cost.py > remaining-cost.csv
"""

import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
RUNS = {
    "QOCGUM_scalar": "keyframe-lean-audit-1250649/QOCGUM-audit",
    "LOFLEA_scalar": "remaining-cost-1250663/LOFLEA-scalar",
    "LOFLEA_two_macro_phases": "remaining-cost-1250665/LOFLEA-phases",
}
SCOPES = {
    "structure_eigensolve": ("ciah_forward_cost", ("ci_solve",)),
    "coupled_expansion_action": ("ciah_expansion_cost", ("hvp_seconds",)),
    "expansion_metric": ("ciah_expansion_cost", ("metric_seconds",)),
    "diagonal_build": ("ciah_diagonal_cost", ("seconds",)),
    "precondition_apply": ("ciah_precondition_cost", ("seconds",)),
    "ah_ritz": ("ciah_ritz_cost", (
        "projection_seconds", "spectrum_image_seconds", "norm_seconds")),
    "keyframe": ("ciah_keyframe_cost", ("seconds",)),
}

writer = csv.writer(sys.stdout)
writer.writerow(("run", "scope", "calls", "wall_seconds"))
for name, relative in RUNS.items():
    events = {}
    phases = {}
    for line in (ROOT / relative / "run.err").read_text().splitlines():
        fields = line.split(",")
        if fields[0] == "NEO_DIAG":
            values = dict(field.split("=", 1) for field in fields[4:] if "=" in field)
            events.setdefault(fields[3], []).append(values)
        elif fields[0] == "XMVB_PHASE" and fields[1] != "stage":
            phases[fields[1]] = (int(fields[2]), float(fields[3]))
    totals = {}
    for scope, (event, keys) in SCOPES.items():
        samples = [sum(float(row.get(key, 0)) for key in keys)
                   for row in events.get(event, [])]
        if scope == "structure_eigensolve":
            samples = [value for value in samples if value > 0]
        totals[scope] = sum(samples)
        writer.writerow((name, scope, len(samples), f"{sum(samples):.9f}"))
    for scope, (calls, seconds) in phases.items():
        writer.writerow((name, "phase_" + scope, calls, f"{seconds:.9f}"))
    if "ci_sigma" in phases:
        lower_bound = max(0, totals["structure_eigensolve"] - phases["ci_sigma"][1])
        writer.writerow((name, "eigensolve_outside_sigma_lower_bound", "",
                         f"{lower_bound:.9f}"))
