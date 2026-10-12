#!/usr/bin/env python3
"""Summarize bounded carrier-connection diagnostics, without reading by hand."""
import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    groups = defaultdict(list)
    endings = {}
    for line in args.log.read_text().splitlines():
        fields = line.split(",")
        if len(fields) < 5 or fields[0] != "NEO_DIAG":
            continue
        event = fields[3]
        if event not in {"ciah_connection", "ciah_micro_end"}:
            continue
        values = dict(item.split("=", 1) for item in fields[4:])
        macro = int(fields[1])
        if event == "ciah_micro_end":
            endings[macro] = values["reason"]
        else:
            groups[macro].append({key: float(value) for key, value in values.items()})
    rows = []
    for macro, records in sorted(groups.items()):
        median = lambda fn: statistics.median(fn(record) for record in records)
        rows.append({
            "macro": macro, "directions": len(records),
            "micro_exit": endings.get(macro, "missing"),
            "median_representation_over_total": median(
                lambda r: r["representation"] / max(r["total"], 1e-300)),
            "median_integral_over_total": median(
                lambda r: r["integral"] / max(r["total"], 1e-300)),
            "median_cross_cosine": median(lambda r: r["cross"] /
                max(r["representation"] * r["integral"], 1e-300)),
            "max_absolute_relative_closure_squared": max(
                abs(r["closure_squared"]) / max(r["carrier_squared"], 1e-300)
                for r in records),
            "median_kinematic_relative_error": median(lambda r:
                r["kinematic_error"] / max(r["representation"], 1e-300)),
            "max_reconstruction_error": max(r["reconstruction"] for r in records),
            "max_metric_relative_residual": max(r["metric_residual"] for r in records),
            "connection_metric_actions": sum(int(r["metric_actions"]) for r in records),
        })
    if not rows:
        raise SystemExit("no ciah_connection events found")
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} macros, {sum(r['directions'] for r in rows)} directions -> {args.output}")


if __name__ == "__main__":
    main()
