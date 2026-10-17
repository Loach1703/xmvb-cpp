#!/usr/bin/env python3
"""Extract opt-in COV stage clocks and work counts from saved Slurm artifacts.

The first-step scope includes initial preparation and the first macro, but not
the next accepted-point diagonal. Nested phase clocks must not be added.
"""

import argparse
from collections import defaultdict
import csv
from pathlib import Path
import sys


def summarize(root):
    totals = defaultdict(lambda: [0, 0.0, 0, 0.0])
    phases, cpu, first_events = {}, {}, defaultdict(lambda: [0, 0.0])
    macro, inside = 0, False

    def add(kind, name, field, value, first):
        row = totals[kind, name, field]
        row[0] += 1
        row[1] += value
        if first:
            row[2] += 1
            row[3] += value

    log = root / "current" / "run.err"
    for line in log.read_text().splitlines():
        fields = line.split(",")
        first = macro == 0 or (macro == 1 and inside)
        if fields[0] == "XMVB_PHASE_EVENT":
            name, event = fields[1:3]
            if name == "ciah_macro" and event == "begin":
                macro += 1
                inside = True
            elif event == "end" and macro == 1 and inside:
                first_events[name][0] += 1
                first_events[name][1] += float(fields[3])
            if name == "ciah_macro" and event == "end":
                inside = False
        elif fields[0] == "XMVB_PHASE" and fields[2].isdigit():
            phases[fields[1]] = int(fields[2]), float(fields[3])
        elif fields[0] == "XMVB_CPU" and fields[1] != "stage":
            cpu[fields[1]] = float(fields[2])
        elif fields[0] == "CIAH_DIAG" and fields[3] in (
            "ciah_diagonal_cost", "ciah_diagonal_sections", "ciah_forward_cost",
            "ciah_structure_cost", "ciah_keyframe_cost", "ciah_keyframe_active_cost",
        ):
            for field in fields[4:]:
                name, value = field.split("=", 1)
                add("scalar", fields[3], name, float(value), first)
        elif fields[0] == "XMVB_TRANSPORT_SOLVE":
            record = dict(field.split("=", 1) for field in fields[1:])
            for name in ("wall", "action_wall", "actions", "iterations", "certifications"):
                add("transport", record["mode"], name, float(record[name]), first)
        elif fields[0] == "XMVB_PRODUCT_ACTION":
            record = dict(field.split("=", 1) for field in fields[1:])
            action = "M" if record["overlap"] == "1" else "H"
            for name in ("rhs", "input", "output", "input_evictions", "output_evictions"):
                add("product", action, name, float(record[name]), first)
        elif fields[0] == "XMVB_TILE":
            record = dict(field.split("=", 1) for field in fields[1:])
            for name in ("prepare", "read", "compute", "write", "source_elements"):
                add("tile", "sigma", name, float(record[name]), first)

    node = (root / "node.txt").read_text().strip()
    for name, (calls, wall) in sorted(phases.items()):
        if calls:
            yield dict(run=root.name, node=node, kind="phase", name=name,
                       field="wall_seconds", calls=calls, total=wall,
                       first_calls=first_events[name][0] if name in first_events else "",
                       first_total=first_events[name][1] if name in first_events else "",
                       cpu_seconds=cpu.get(name, ""))
    for (kind, name, field), (calls, total, first_calls, first_total) in sorted(totals.items()):
        yield dict(run=root.name, node=node, kind=kind, name=name, field=field,
                   calls=calls, total=total, first_calls=first_calls,
                   first_total=first_total, cpu_seconds="")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path, nargs="+")
    args = parser.parse_args()
    rows = [row for root in args.runs for row in summarize(root)]
    if rows:
        writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
