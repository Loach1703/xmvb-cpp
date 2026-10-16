#!/usr/bin/env python3
"""Summarize saved one-step GitHub comparisons, including failed/unfinished runs."""

import argparse
import csv
from pathlib import Path
import re
import sys


def field(text, pattern, convert=str):
    found = re.search(pattern, text, re.MULTILINE)
    return convert(found[1]) if found else None


def elapsed_seconds(text):
    value = field(text, r"Elapsed \(wall clock\).*?: ([0-9:.]+)")
    if value is None:
        return None
    total = 0.0
    for part in value.split(":"):
        total = 60 * total + float(part)
    return total


def summarize(root):
    space = (root / "space.txt").read_text().strip()
    node = (root / "node.txt").read_text().strip()
    for variant in ("baseline", "current"):
        run = root / variant
        output = (run / "output.out").read_text() if (run / "output.out").exists() else ""
        resource = (run / "resource.time").read_text() if (run / "resource.time").exists() else ""
        trace = []
        if (run / "neo.tsv").exists():
            with (run / "neo.tsv").open() as stream:
                trace = list(csv.DictReader(stream, delimiter="\t"))
        yield dict(
            job=root.name, node=node, space=space, variant=variant,
            actual_threads=8,
            exit_code=(run / "run.exit").read_text().strip() if (run / "run.exit").exists() else None,
            termination=field(output, r"^Termination reason\s*:\s*(\S+)"),
            structures=field(output, r"Number of structures:\s*(\d+)", int),
            recorded_accepted_steps=len(trace),
            orbital_hvp=sum(int(row["orbital_hvp_actions"]) for row in trace) if trace else None,
            structure_response=sum(int(row["structure_response_actions"]) for row in trace) if trace else None,
            initial_Eh=field(output, r"^Initial total energy\s*:\s*([-\d.]+)", float),
            after_step_Eh=field(output, r"^Final total energy\s*:\s*([-\d.]+)", float),
            scf_s=field(output, r"^SCF iteration wall time\s*:\s*([\d.]+)", float),
            application_s=field(output, r"^End-to-end wall time\s*:\s*([\d.]+)", float),
            process_s=elapsed_seconds(resource),
            cpu_s=(field(resource, r"User time \(seconds\):\s*([\d.]+)", float) +
                   field(resource, r"System time \(seconds\):\s*([\d.]+)", float))
                if "User time (seconds)" in resource else None,
            rss_kib=field(resource, r"Maximum resident set size \(kbytes\):\s*(\d+)", int),
            cpu_percent=field(resource, r"Percent of CPU this job got:\s*(\d+)%", int),
            last_ci_actions=field(output, r"^Final Davidson H/S block actions\s*:\s*(\d+)", int),
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="+", type=Path)
    args = parser.parse_args()
    rows = [row for run in args.runs for row in summarize(run)]
    writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
