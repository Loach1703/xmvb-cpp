#!/usr/bin/env python3
"""Summarize existing native C12 artifacts without launching calculations."""
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "tests/mcscf"))
from summarize_scf_step import elapsed_seconds, field

RUNS = [
    ("saved_scatter", "c12-scatter-address-20261017", 1264636),
    ("saved_halo", "c12-halo-lifetime-20261017", 1264676),
    ("saved_halo_repeat", "c12-halo-lifetime-20261017", 1264685),
    ("three_slab_unwired", HERE.name, 1264695),
    ("ordered_two_slab_unwired", HERE.name, 1264696),
    ("joint_physical", HERE.name, 1264713),
]


def row(label, directory, job):
    run = HERE.parent / directory / f"run-{job}" / "native"
    output = (run / "output.out").read_text()
    resource = (run / "resource.time").read_text()
    with (run / "ciah.tsv").open() as stream:
        trace = list(csv.DictReader(stream, delimiter="\t"))
    return dict(
        variant=label, job=job,
        steps=len(trace),
        hvp=sum(int(r["orbital_hvp_actions"]) for r in trace),
        structure_responses=sum(int(r["structure_response_actions"]) for r in trace),
        energy_Eh=field(output, r"^Final total energy\s*:\s*([-\d.]+)"),
        scf_s=field(output, r"^SCF.*wall.*:\s*([\d.]+)"),
        application_s=field(output, r"^End-to-end wall time\s*:\s*([\d.]+)"),
        process_s=elapsed_seconds(resource),
        rss_MiB=field(resource, r"Maximum resident set size \(kbytes\):\s*(\d+)", int)/1024,
        termination=field(output, r"^Termination reason\s*:\s*(\S+)"))


if __name__ == "__main__":
    rows = [row(*run) for run in RUNS]
    with (HERE / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
