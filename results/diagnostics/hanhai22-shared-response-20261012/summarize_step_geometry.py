#!/usr/bin/env python3
"""Extract actual orbital-step geometry; no chemistry or thresholds are changed."""
import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path


def quantile(values, fraction):
    values = sorted(values)
    position = fraction * (len(values) - 1)
    lower = int(position)
    upper = min(lower + 1, len(values) - 1)
    return values[lower] + (position - lower) * (values[upper] - values[lower])


def read_geometry(directory):
    pattern = re.compile(r"^NEO_DIAG,(\d+),(\d+),orbital_step_geometry,(.*)$")
    rows = []
    batch = 0
    with (directory / "run.err").open() as stream:
        for line_number, line in enumerate(stream, 1):
            match = pattern.match(line.strip())
            if not match:
                continue
            fields = dict(item.split("=", 1) for item in match[3].split(","))
            if int(fields["orbital"]) == 1:
                batch += 1
            row = dict(case=directory.name, macro=int(match[1]), attempt=int(match[2]),
                       batch=batch, stage=fields["stage"], micro=int(fields["micro"]),
                       orbital=int(fields["orbital"]), role=fields["role"],
                       coefficients=int(fields["coefficients"]), source_line=line_number)
            for field in ("alpha", "theta", "step_norm", "actual_angle", "norm_ratio"):
                row[field] = float(fields[field])
            row["radial_fraction"] = (abs(row["alpha"]) / row["step_norm"]
                                      if row["step_norm"] else 0.0)
            row["tangent_add_angle"] = math.atan(row["theta"])
            row["radial_angle_change"] = row["actual_angle"] - row["tangent_add_angle"]
            row["actual_minus_exponential_angle"] = row["actual_angle"] - row["theta"]
            rows.append(row)
    return rows


def summarize(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[row["case"], row["stage"], row["role"]].append(row)
    result = []
    for (case, stage, role), group in sorted(groups.items()):
        radial_squared = sum(row["alpha"] ** 2 for row in group)
        total_squared = sum(row["step_norm"] ** 2 for row in group)
        theta = [row["theta"] for row in group]
        angle = [row["actual_angle"] for row in group]
        # Weighted fraction avoids giving tiny displacements undue importance.
        result.append(dict(
            case=case, stage=stage, role=role, samples=len(group),
            batches=len({row["batch"] for row in group}),
            macros=len({row["macro"] for row in group}),
            alpha_min=min(row["alpha"] for row in group),
            alpha_max=max(row["alpha"] for row in group),
            alpha_abs_max=max(abs(row["alpha"]) for row in group),
            radial_l2_fraction=math.sqrt(radial_squared / total_squared) if total_squared else 0.0,
            theta_median=quantile(theta, .5), theta_p95=quantile(theta, .95), theta_max=max(theta),
            actual_angle_p95=quantile(angle, .95), actual_angle_max=max(angle),
            max_abs_radial_angle_change=max(abs(row["radial_angle_change"]) for row in group),
            max_abs_actual_minus_exponential_angle=max(
                abs(row["actual_minus_exponential_angle"]) for row in group),
            norm_ratio_min=min(row["norm_ratio"] for row in group),
            norm_ratio_max=max(row["norm_ratio"] for row in group)))
    return result


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directories", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    rows = [row for directory in arguments.directories for row in read_geometry(directory)]
    if not rows:
        raise SystemExit("No orbital_step_geometry records found")
    arguments.output.mkdir(parents=True, exist_ok=True)
    summary = summarize(rows)
    write_csv(arguments.output / "orbital_steps.csv", rows)
    write_csv(arguments.output / "accepted_orbital_steps.csv",
              [row for row in rows if row["stage"] == "macro_accepted"])
    write_csv(arguments.output / "geometry_summary.csv", summary)
    (arguments.output / "geometry_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for row in summary:
        print(f'{row["case"]:8s} {row["stage"]:15s} {row["role"]:8s} '
              f'n={row["samples"]:5d} max|alpha|={row["alpha_abs_max"]:.6g} '
              f'theta_p95={row["theta_p95"]:.6g} theta_max={row["theta_max"]:.6g} '
              f'radial_l2={row["radial_l2_fraction"]:.4%} '
              f'max|angle-theta|={row["max_abs_actual_minus_exponential_angle"]:.6g}')


if __name__ == "__main__":
    main()
