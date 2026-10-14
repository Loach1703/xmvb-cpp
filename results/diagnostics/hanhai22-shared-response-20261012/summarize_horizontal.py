"""Summarize saved S-horizontal regressions; never launches calculations."""
import csv
import gzip
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent


def measure(folder):
    text = (folder / "output.out").read_text()
    resource = (folder / "resource.time").read_text()
    rows = list(csv.DictReader((folder / "neo.tsv").open(), delimiter="\t"))
    wall = re.search(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)", resource).group(1)
    seconds = sum(float(part) * 60**i for i, part in enumerate(wall.split(":")[::-1]))
    return dict(iterations=len(rows), hvp=sum(int(r["orbital_hvp_actions"]) for r in rows),
                energy=float(re.findall(r"TOTAL ENERGY\s*:\s*([-\d.]+)", text)[-1]),
                wall_seconds=seconds,
                rss_mib=int(re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", resource).group(1))/1024,
                exit_code=int((folder / "run.exit").read_text()),
                converged="VBSCF converged in" in text)


def trace(folder):
    result = []
    for line in gzip.open(folder / "run.err.gz", "rt"):
        fields = line.strip().split(",")
        if len(fields) >= 4 and fields[0] == "NEO_DIAG":
            values = dict(field.split("=", 1) for field in fields[4:] if "=" in field)
            result.append(dict(macro=int(fields[1]), event=fields[3], **values))
    return result


def main():
    comparison = []
    for case in ("F2-OEO", "FeCl2", "MnF2", "241"):
        before = measure(ROOT / "transport-regression-1254179" / case)
        after = measure(ROOT / "horizontal-regression-1254764" / case)
        assert after["converged"] and after["exit_code"] == 0
        row = dict(case=case, delta_energy=after["energy"]-before["energy"])
        row.update({"before_"+k: v for k, v in before.items()})
        row.update({"after_"+k: v for k, v in after.items()})
        comparison.append(row)
    report = dict(molecular_regression=comparison,
                  singlet_short=measure(ROOT / "horizontal-regression-1254764/singlet"))
    report["singlet_short"]["purpose"] = "Two-macro derivative diagnostic, not convergence or performance"
    records = trace(ROOT / "horizontal-regression-1254764/singlet")
    report["singlet_derivative_audits"] = [r for r in records if r["event"] in
        ("ciah_transported_hvp_audit", "ciah_transported_gradient_audit")]
    report["singlet_keyframes"] = [r for r in records if r["event"] in
        ("ciah_keyframe", "ciah_keyframe_blocks", "ciah_macro_trial")]
    report["raw_active_geometry"] = []
    for macro in (1, 2):
        for stage in ("keyframe", "macro_accepted"):
            selected = [r for r in records if r["macro"] == macro and
                r["event"] == "orbital_step_geometry" and r["stage"] == stage and r["role"] == "active"]
            report["raw_active_geometry"].append(dict(macro=macro, stage=stage,
                max_abs_alpha=max(abs(float(r["alpha"])) for r in selected),
                max_theta=max(float(r["theta"]) for r in selected)))
    for label, folder in (("before_geometry", "representative-analysis-1254511"),
                          ("after_geometry", "representative-analysis-1254782")):
        data = json.loads((ROOT / folder / "representative-analysis.json").read_text())
        selected = next(v for k, v in data["candidate_geometry"].items() if k.startswith("macro-002-keyframe"))
        report[label] = dict(
            active_angle=selected["geometry"]["active"]["max_angle_degrees"],
            physical_norm=selected["independent_physical_tangent_norm"],
            parallel_norm=selected["active_frame_tangent"]["parallel_frobenius"],
            orientation=selected["active_finite_frame_rotation"]["determinant"],
            isometry_difference=selected["raw_vs_same_tangent_isometry"]["active"]["max_angle_degrees"],
            native_pyscf_macro2=data["between"]["neo2_pyscf2"]["active"]["max_angle_degrees"])
    output = ROOT / "horizontal-regression-1254764/summary"
    output.mkdir(exist_ok=True)
    (output / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    with (output / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=comparison[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(comparison)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
