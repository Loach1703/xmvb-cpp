"""Summarize scalar CIAH accounting; no manual trajectory transcription."""
import csv
import json
import pathlib
import sys

run = pathlib.Path(sys.argv[1])
events = {}
ritz_blocks = {}
last_norm = None
for line in (run / "run.err").read_text().splitlines():
    if not line.startswith("NEO_DIAG,"):
        continue
    _, macro, _, event, payload = line.split(",", 4)
    values = dict(item.split("=", 1) for item in payload.split(",") if "=" in item)
    events.setdefault(event, []).append((int(macro), values))
    if event == "ciah_coupled_norm":
        last_norm = (int(macro), values)
    elif event == "ciah_ritz_cost" and last_norm and last_norm[0] == int(macro):
        ritz_blocks.setdefault(int(macro), []).append(last_norm[1])
norms = [v for _, v in events.get("ciah_structure_norm", [])]
ritz = [v for _, v in events.get("ciah_ritz_cost", [])]
expansions = [v for _, v in events.get("ciah_expansion_cost", [])]
preconditions = [v for _, v in events.get("ciah_precondition_cost", [])]
forward = [v for _, v in events.get("ciah_forward_cost", [])]
diagonals = [v for _, v in events.get("ciah_diagonal_cost", [])]
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
    "precondition_seconds": sum(float(v["seconds"]) for v in preconditions),
    "diagonal_builds": len(diagonals),
    "diagonal_build_seconds": sum(float(v["seconds"]) for v in diagonals),
    "forward_evaluations": len(forward),
    "forward_seconds": {key: sum(float(v[key]) for v in forward) for key in
                        ("orbital", "ao_effective", "active_one", "active_two",
                         "structure_setup", "ci_solve")},
    "expansion_seconds": {key: sum(float(v[key]) for v in expansions) for key in
                          ("orthogonalize_seconds", "metric_seconds", "hvp_seconds")},
    "max_relative_galerkin_residual": max(
        (float(v["galerkin_residual"]) /
         max(float(v["residual_euclidean"]), 1.0e-300) for v in ritz), default=0),
    "ritz_seconds": {key: sum(float(v[key]) for v in ritz) for key in
                     ("projection_seconds", "spectrum_image_seconds", "norm_seconds")},
    "keyframe_direction_audits": [
        {"macro": macro, "delta": float(v["delta"]),
         "derivative_error": float(v["derivative_error"]),
         "relative_curvature_change": float(v["curvature_change"]) /
             max(float(v["origin_curvature"]), 1.0e-300)}
        for macro, v in events.get("ciah_keyframe_direction_audit", [])],
    "keyframe_secants": [
        {"macro": macro, "mismatch": float(v["mismatch"]),
         "captured_squared_fraction": (float(v["captured"]) /
             max(float(v["mismatch"]), 1.0e-300))**2,
         "secant_cosine": float(v["secant_cosine"])}
        for macro, v in events.get("ciah_keyframe_secant", [])],
    "schur_spectra": [dict(macro=macro, **v)
                      for macro, v in events.get("ciah_schur_spectrum", [])],
    "metric_decision_disagreements": sum(
        v["diagonal_keep"] != v["full_keep"] or
        v["diagonal_trigger"] != v["full_trigger"]
        for _, v in events.get("ciah_keyframe_metric_audit", [])),
    "tail": [],
}
for macro, end in events.get("ciah_micro_end", [])[-8:]:
    micro = [v for k, v in events.get("ciah_micro", []) if k == macro]
    keyframes = [v for k, v in events.get("ciah_keyframe", []) if k == macro]
    additions = [v for k, v in events.get("ciah_expansion_cost", []) if k == macro]
    blocks = [v for k, v in events.get("ciah_block_curvature", []) if k == macro]
    corrections = [v for k, v in events.get("ciah_keyframe_blocks", []) if k == macro]
    summary["tail"].append({
        "macro": macro, "stop": end["reason"],
        "micro": len(micro), "source": float(end["source"]),
        "minimum_predicted": min(float(v["gradient"]) for v in micro),
        "last_predicted": float(end["predicted"]),
        "clipped": sum(float(v["step_scale"]) < 1 for v in micro),
        "last_residual": float(micro[-1]["residual"]) if micro else None,
        "last_ah_eigenvalue": float(micro[-1]["omega"]) if micro else None,
        "last_ritz_norm_blocks": ritz_blocks.get(macro, [None])[-1],
        "keyframes": len(keyframes),
        "minimum_new_direction_fraction": min(
            (float(v["remaining_fraction"]) for v in additions), default=None),
        "max_keyframe_relative_mismatch": max(
            (float(v["mismatch"]) / max(float(v["predicted"]), 1.0e-300)
             for v in keyframes), default=0),
        "max_coupling_fraction": max(
            (abs(float(v["coupling_forward"])+float(v["coupling_backward"])) /
             max(abs(float(v["orbital"]))+abs(float(v["structure"])),1.0e-300)
             for v in blocks), default=None),
        "min_positive_diagonal_curvature_retained": min(
            ((float(v["orbital"])+float(v["structure"])+
              float(v["coupling_forward"])+float(v["coupling_backward"])) /
             (float(v["orbital"])+float(v["structure"]))
             for v in blocks
             if float(v["orbital"]) > 0 and float(v["structure"]) > 0), default=None),
        "keyframe_mismatch_orbital_fractions": [
            float(v["mismatch_orbital"])**2 /
            max(float(v["mismatch_orbital"])**2+float(v["mismatch_structure_squared"]),1.0e-300)
            for v in corrections],
    })
print(json.dumps(summary, indent=2))
