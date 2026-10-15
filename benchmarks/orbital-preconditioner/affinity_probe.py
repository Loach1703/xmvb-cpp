#!/usr/bin/env python3
"""Record the allocated CPU topology and launch on a nested physical-core set."""
import argparse
import json
import os
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--threads", type=int, required=True)
parser.add_argument("--nested", action="store_true")
parser.add_argument("--compact", action="store_true")
parser.add_argument("--record", required=True)
parser.add_argument("command", nargs=argparse.REMAINDER)
args = parser.parse_args()
allowed = sorted(os.sched_getaffinity(0))
topology = []
for cpu in allowed:
    base = Path(f"/sys/devices/system/cpu/cpu{cpu}")
    numa = [p.name for p in base.glob("node*")]
    topology.append(dict(cpu=cpu, numa=numa,
                         core=(base / "topology/core_id").read_text().strip(),
                         socket=(base / "topology/physical_package_id").read_text().strip(),
                         l3=(base / "cache/index3/shared_cpu_list").read_text().strip()))
if args.compact:
    # Stay inside the scheduler's allowed mask. Prefer the largest NUMA domain,
    # then complete L3 domains; never request CPUs belonging to another job.
    groups = {}
    for item in topology:
        numa = tuple(item["numa"])
        groups.setdefault(numa, {}).setdefault(item["l3"], []).append(item["cpu"])
    ordered = []
    for _, caches in sorted(groups.items(), key=lambda item: (-sum(map(len, item[1].values())), item[0])):
        for _, cpus in sorted(caches.items(), key=lambda item: (-len(item[1]), item[1][0])):
            ordered.extend(cpus)
    selected = ordered[:args.threads]
else:
    selected = allowed[:args.threads] if args.nested else allowed
if len(selected) < args.threads:
    raise RuntimeError("allocation has fewer CPUs than requested threads")
Path(args.record).write_text(json.dumps(dict(allowed=allowed, selected=selected,
                                           topology=topology), indent=2) + "\n")
os.sched_setaffinity(0, selected)
command = args.command
if command and command[0] == "--":
    command = command[1:]
os.execvp(command[0], command)
