#!/usr/bin/env python3
"""Extract phase CPU accounting and worker balance without summing nested scopes."""
import argparse
import csv
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("root", type=Path)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
rows = []
for directory in sorted(args.root.iterdir()):
    log = directory / "run.err"
    if not log.is_file():
        continue
    resources = directory / "resource.time"
    if not resources.is_file() or "Exit status: 0" not in resources.read_text():
        continue
    phases, cpu, balance = {}, {}, {}
    for line in log.read_text().splitlines():
        fields = line.split(",")
        if len(fields) < 3 or fields[1] == "stage":
            continue
        if fields[0] == "XMVB_PHASE":
            phases[fields[1]] = (int(fields[2]), float(fields[3]))
        elif fields[0] == "XMVB_CPU":
            cpu[fields[1]] = float(fields[2])
        elif fields[0] == "XMVB_BALANCE":
            balance[fields[1]] = list(map(float, fields[2:]))
    threads = int(directory.name.rsplit("-", 1)[1])
    for stage, (calls, wall) in phases.items():
        if not calls:
            continue
        row = dict(variant=directory.name, threads=threads, stage=stage,
                   calls=calls, wall_seconds=wall, cpu_seconds=cpu.get(stage, ""),
                   cpu_wall_ratio=cpu[stage] / wall if wall and stage in cpu else "",
                   loop_occupancy="", loop_min_max_ratio="", start_skew_seconds="")
        if stage in balance:
            _, minimum, maximum, total, capacity, skew, _ = balance[stage]
            row.update(loop_occupancy=total / capacity if capacity else "",
                       loop_min_max_ratio=minimum / maximum if maximum else "",
                       start_skew_seconds=skew)
        rows.append(row)
if not rows:
    raise RuntimeError("no completed phase logs found")
args.output.parent.mkdir(parents=True, exist_ok=True)
with args.output.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(f"Saved {len(rows)} phase records to {args.output}")
print("Scopes overlap: do not sum their times. CPU includes spinning; short-window accounting is approximate.")
print("Use whole-process GNU-time averages, not phase CPU/wall ratios, for utilization.")
