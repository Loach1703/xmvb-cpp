#!/usr/bin/env python3
"""Compare saved 240 FULL runs and nested phase records; never runs binaries."""

import argparse
import csv
import hashlib
from pathlib import Path
import re


def scalar(text, pattern):
    match = re.search(pattern, text, re.MULTILINE)
    if not match:
        raise ValueError(f"missing field: {pattern}")
    return match[1]


def phases(root):
    wall, cpu = {}, {}
    for line in (root / "phases/full/run.err").read_text().splitlines():
        fields = line.split(",")
        if fields[0] == "XMVB_PHASE" and fields[1] != "stage":
            wall[fields[1]] = (int(fields[2]), float(fields[3]))
        elif fields[0] == "XMVB_CPU" and fields[1] != "stage":
            cpu[fields[1]] = float(fields[2])
    return wall, cpu


def native(root):
    case = root / "native/full"
    text = (case / "output.out").read_text()
    resource = (case / "resource.time").read_text()
    trace = list(csv.DictReader((case / "neo.tsv").open(), delimiter="\t"))
    seconds = 0.0
    for field in scalar(resource, r"Elapsed \(wall clock\).*?: ([0-9:.]+)").split(":"):
        seconds = 60 * seconds + float(field)
    user = float(scalar(resource, r"User time \(seconds\): ([0-9.]+)"))
    system = float(scalar(resource, r"System time \(seconds\): ([0-9.]+)"))
    return dict(
        run=root.name, mode="native", space="full",
        node=(root / "node.txt").read_text().strip(),
        threads=int(scalar((root / "environment.txt").read_text(), r"^OMP_NUM_THREADS=(\d+)")),
        steps=int(scalar(text, r"VBSCF converged in\s+(\d+) iterations")),
        hvp=sum(int(row["orbital_hvp_actions"]) for row in trace),
        energy_Eh=float(scalar(text, r"Final total energy\s+:\s+([-0-9.]+)")),
        scf_s=float(scalar(text, r"SCF iteration wall time\s+:\s+([0-9.]+)")),
        app_s=float(scalar(text, r"End-to-end wall time\s+:\s+([0-9.]+)")),
        process_s=seconds, rss_mib=float(scalar(resource, r"Maximum resident set size \(kbytes\):\s*(\d+)"))/1024,
        process_cpu_s=user+system, mean_process_cpu_cores=(user+system)/seconds)


def save(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("historical", type=Path)
    parser.add_argument("current", type=Path)
    parser.add_argument("historical_summary", type=Path)
    parser.add_argument("output_prefix", type=Path)
    args = parser.parse_args()
    old = next(row for row in csv.DictReader(args.historical_summary.open())
               if row["job"] == args.historical.name.rsplit("-", 1)[1]
               and row["mode"] == "native" and row["space"] == "full")
    new = native(args.current)
    old_row = {key: old.get(key, "") for key in new}
    old_threads = int(scalar((args.historical / "environment.txt").read_text(), r"^OMP_NUM_THREADS=(\d+)"))
    old_row.update(run=args.historical.name, threads=old_threads)
    old_time = (args.historical / "native/full/resource.time").read_text()
    old_row["process_cpu_s"] = sum(float(scalar(old_time, key + r": ([0-9.]+)"))
        for key in (r"User time \(seconds\)", r"System time \(seconds\)"))
    old_row["mean_process_cpu_cores"] = old_row["process_cpu_s"] / float(old_row["process_s"])
    save(args.output_prefix.with_suffix(".csv"), [old_row, new])
    old_wall, old_cpu = phases(args.historical)
    new_wall, new_cpu = phases(args.current)
    rows = []
    for stage in sorted(old_wall.keys() | new_wall.keys()):
        oc, ow = old_wall.get(stage, (0, 0.0))
        nc, nw = new_wall.get(stage, (0, 0.0))
        rows.append(dict(stage=stage, historical_calls=oc, current_calls=nc,
                         historical_wall_s=ow, current_wall_s=nw, delta_wall_s=nw-ow,
                         historical_cpu_s=old_cpu.get(stage, ""), current_cpu_s=new_cpu.get(stage, "")))
    save(args.output_prefix.with_name(args.output_prefix.name + "_phases.csv"), rows)
    for filename in ("input.xmi", "cc-pvdz.gbs"):
        digests = [hashlib.sha256((root/"native/full"/filename).read_bytes()).hexdigest()
                   for root in (args.historical, args.current)]
        print(filename, "identical" if digests[0] == digests[1] else "DIFFERENT", *digests)
    print("Native current:", new)
    print("Largest positive nested phase differences (not additive):")
    for row in sorted(rows, key=lambda row: row["delta_wall_s"], reverse=True)[:15]:
        print(row)


if __name__ == "__main__":
    main()
