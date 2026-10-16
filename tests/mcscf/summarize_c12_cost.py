#!/usr/bin/env python3
"""Parse saved C12 phase probes and scalar CIAH events; never run a workload."""

import argparse
import collections
import csv
import json
from pathlib import Path
import statistics

from summarize_scf_step import elapsed_seconds, field


def values(parts):
    return {key: float(value) for key, value in
            (part.split("=", 1) for part in parts if "=" in part)}


def summarize(run):
    case = run / "loflea-phases"
    output = (case / "output.out").read_text()
    resource = (case / "resource.time").read_text()
    phases, cpus, tiles = {}, {}, []
    events = collections.defaultdict(list)
    for line in (case / "run.err").read_text().splitlines():
        parts = line.split(",")
        if parts[0] == "XMVB_PHASE" and parts[2].isdigit():
            phases[parts[1]] = dict(calls=int(parts[2]), wall_s=float(parts[3]))
        elif parts[0] == "XMVB_CPU" and parts[1] != "stage":
            cpus[parts[1]] = float(parts[2])
        elif parts[0] == "XMVB_TILE":
            tiles.append(values(parts[1:]))
        elif parts[0] == "CIAH_DIAG":
            if parts[3] in ("ciah_structure_cost", "ciah_transport_magnitude"):
                events[parts[3]].append(values(parts[4:]))
    for name, data in phases.items():
        data["process_cpu_s"] = cpus.get(name)
    with (case / "ciah.tsv").open() as stream:
        trace = list(csv.DictReader(stream, delimiter="\t"))
    event_stats = {}
    for event, records in events.items():
        event_stats[event] = dict(calls=len(records), fields={})
        for key in records[0]:
            samples = [row[key] for row in records]
            event_stats[event]["fields"][key] = dict(
                sum=sum(samples), min=min(samples),
                median=statistics.median(samples), max=max(samples))
    tile_totals = dict(sweeps=len(tiles), **{
        key: sum(row[key] for row in tiles) for key in
        ("prepare", "read", "compute", "write", "source_elements", "batches")})
    # Older saved runs did not record these scopes: keep them missing, not zero.
    for key in ("gram", "gram_tiles"):
        measured = [row[key] for row in tiles if key in row]
        if measured:
            tile_totals[key] = sum(measured)
    return dict(
        run=run.name, node=(run / "node.txt").read_text().strip(),
        note="Nested phase totals overlap; they are not additive wall-time shares.",
        initial_Eh=field(output, r"^Initial total energy\s*:\s*([-\d.]+)", float),
        final_Eh=field(output, r"^Final total energy\s*:\s*([-\d.]+)", float),
        termination=field(output, r"^Termination reason\s*:\s*(\S+)"),
        steps=len(trace), hvp=sum(int(row["orbital_hvp_actions"]) for row in trace),
        application_s=field(output, r"^End-to-end wall time\s*:\s*([\d.]+)", float),
        process_s=elapsed_seconds(resource),
        rss_kib=field(resource, r"Maximum resident set size \(kbytes\):\s*(\d+)", int),
        phases=phases,
        tile_totals=tile_totals,
        events=event_stats)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps([summarize(run) for run in args.runs], indent=2) + "\n"
    if args.output:
        args.output.write_text(result)
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
