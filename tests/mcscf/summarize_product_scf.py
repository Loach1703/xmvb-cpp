#!/usr/bin/env python3
"""Parse current converged 240 runs and compare saved equal-thread FULL controls."""

import argparse
import csv
import json
from pathlib import Path
import re


def scalar(text, pattern):
    match = re.search(pattern, text, re.MULTILINE)
    if not match:
        raise ValueError(f"missing {pattern}")
    return match[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("saved", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with (args.saved / "selected_transport_240.csv").open() as stream:
        reference = {row["space"]:row for row in csv.DictReader(stream)}
    with (args.saved / "full32_comparison.csv").open() as stream:
        controls = list(csv.DictReader(stream))
    rows=[]
    topology=args.run/"topology.json"
    sockets=len(json.loads(topology.read_text())["sockets"]) if topology.exists() else None
    for space in ("cov","ion01","ion02","ion03","full"):
        run=args.run / "native" / space
        text=(run / "output.out").read_text()
        resource=(run / "resource.time").read_text()
        with (run / "neo.tsv").open() as stream:
            trace=list(csv.DictReader(stream,delimiter="\t"))
        energy=float(scalar(text,r"Final total energy\s*:\s*([-\d.]+)"))
        elapsed=0.0
        for field in scalar(resource,r"Elapsed \(wall clock\).*?: ([0-9:.]+)").split(":"):
            elapsed=60*elapsed+float(field)
        row=dict(run=args.run.name,node=(args.run / "node.txt").read_text().strip(),
                 threads=int(scalar((args.run / "environment.txt").read_text(),r"^OMP_NUM_THREADS=(\d+)")),
                 sockets=sockets,
                 space=space,structures=int(reference[space]["structures"]),
                 steps=int(scalar(text,r"VBSCF converged in\s+(\d+) iterations")),
                 hvp=sum(int(item["orbital_hvp_actions"]) for item in trace),
                 initial_Eh=float(scalar(text,r"Initial total energy\s*:\s*([-\d.]+)")),
                 final_Eh=energy,energy_delta_from_saved_Eh=energy-float(reference[space]["final_energy_Eh"]),
                 scf_s=float(scalar(text,r"SCF iteration wall time\s*:\s*([\d.]+)")),
                 application_s=float(scalar(text,r"End-to-end wall time\s*:\s*([\d.]+)")),
                 process_s=elapsed,rss_mib=int(scalar(resource,r"Maximum resident set size \(kbytes\):\s*(\d+)"))/1024,
                 saved_32core_scf_s=None,saved_32core_rss_mib=None,
                 parent_32core_scf_s=None,parent_32core_rss_mib=None)
        if space=="full":
            row.update(saved_32core_scf_s=float(controls[0]["scf_s"]),
                       saved_32core_rss_mib=float(controls[0]["rss_mib"]),
                       parent_32core_scf_s=float(controls[1]["scf_s"]),
                       parent_32core_rss_mib=float(controls[1]["rss_mib"]))
        rows.append(row)
    with args.output.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    for row in rows:
        print(row)


if __name__=="__main__":
    main()
