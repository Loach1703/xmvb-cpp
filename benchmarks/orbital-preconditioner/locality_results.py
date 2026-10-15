#!/usr/bin/env python3
"""Extract converged-run resources, phase times and placement from saved logs."""
import argparse
import csv
import json
import re
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("root", type=Path)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
rows = []
for directory in sorted(args.root.iterdir()):
    if not (directory / "resource.time").is_file():
        continue
    output = (directory / "output.out").read_text()
    resource = (directory / "resource.time").read_text()
    stderr = (directory / "run.err").read_text()
    iterations = re.search(r"VBSCF converged in\s+(\d+) iterations", output)
    if not iterations or "Exit status: 0" not in resource:
        raise RuntimeError(f"incomplete/unconverged run: {directory}")
    def value(pattern, text):
        match = re.search(pattern, text)
        if not match:
            raise RuntimeError(f"missing field {pattern}: {directory}")
        return float(match.group(1))
    elapsed = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)", resource)
    seconds = 0.0
    for component in elapsed.group(1).split(":"):
        seconds = 60.0 * seconds + float(component)
    cpu = value(r"User time \(seconds\):\s*([\d.]+)", resource) + value(r"System time \(seconds\):\s*([\d.]+)", resource)
    affinity = json.loads((directory / "affinity.json").read_text())
    selected = set(affinity["selected"])
    topology = [item for item in affinity["topology"] if item["cpu"] in selected]
    phases = {}
    calls = {}
    for line in stderr.splitlines():
        if line.startswith("XMVB_PHASE,") and not line.startswith("XMVB_PHASE,stage,"):
            _, phase, count, wall = line.split(",")
            phases[phase] = float(wall)
            calls[phase] = int(count)
    rows.append(dict(variant=directory.name, threads=int(directory.name.split("-")[-1]),
                     steps=int(iterations.group(1)), hvps=calls["hvp_apply"],
                     energy=value(r"Final total energy\s*:\s*([-\d.]+)", output),
                     wall=value(r"End-to-end wall time\s*:\s*([\d.]+)", output),
                     elapsed=seconds, rss_mib=value(r"Maximum resident set size \(kbytes\):\s*(\d+)", resource) / 1024,
                     average_cores=cpu / seconds,
                     numa_domains=len({tuple(item["numa"]) for item in topology}),
                     l3_domains=len({item["l3"] for item in topology}),
                     cpus=" ".join(map(str, affinity["selected"])),
                     **{phase: phases[phase] for phase in ["ao_pair_action", "ao_h1e_fused", "ao_h1e_transpose", "selected_direction", "metric_solve"]}))
if not rows:
    raise RuntimeError("no converged native runs found")
ones = {row["variant"].split("-")[0]: row for row in rows if row["threads"] == 1}
fallback = next(iter(ones.values()), None)
for row in rows:
    one = ones.get(row["variant"].split("-")[0], fallback)
    row["speedup"] = one["wall"] / row["wall"] if one else ""
    row["parallel_efficiency"] = row["speedup"] / row["threads"] if one else ""
args.output.parent.mkdir(parents=True, exist_ok=True)
with args.output.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
for row in rows:
    print(json.dumps(row))
