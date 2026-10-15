#!/usr/bin/env python3
"""Prepare and summarize real selected-structure 240 SCF regressions."""

import argparse
import csv
from pathlib import Path
import re
import shutil


def prepare(bundle, artifact, system):
    text = (bundle / ("240_vb.xmi" if system == "240" else "csd_LOFLEA_vb.xmi")).read_text()
    # 240 has 24 inactive pairs and eight active electrons (orbitals 25--32).
    # These are three actual COV structures, not a completed variational space.
    # List bond endpoints separately: '-' means an integer range in $STR.
    structures = (
        "1:24 28 29 27 30 26 31 25 32",
        "1:24 27 28 29 30 26 31 25 32",
        "1:24 27 28 26 29 30 31 25 32",
    )
    if system == "LOFLEA":
        # Twelve active electrons in orbitals 27--38, after 26 inactive pairs.
        structures = (
            "1:26 27 28 29 30 31 32 33 34 35 36 37 38",
            "1:26 27 30 28 29 31 32 33 34 35 36 37 38",
            "1:26 27 28 29 32 30 31 33 34 35 36 37 38",
        )
    for count in (1, 3):
        deck = re.sub(r"^str=.*\n", "", text, flags=re.I | re.M)
        deck = re.sub(r"^iscf=.*$", "iscf=9", deck, flags=re.I | re.M)
        deck += "\n$STR\n" + "\n".join(structures[:count]) + "\n$END\n"
        for solver in ("davidson", "dense"):
            run = artifact / f"selected{count}-{solver}"
            run.mkdir(parents=True, exist_ok=True)
            (run / "input.xmi").write_text(deck)
            shutil.copyfile(bundle / "cc-pvdz.gbs", run / "cc-pvdz.gbs")


def summarize(artifact, output):
    rows = []
    for count in (1, 3):
        for solver in ("davidson", "dense"):
            run = artifact / f"selected{count}-{solver}"
            deck = (run / "input.xmi").read_text()
            row = dict(
                run=artifact.name, structures=count, eigensolver=solver,
                active_orbitals=int(re.search(r"^nao=(\d+)", deck, re.I | re.M)[1]),
                status="not_started", termination_reason=None, iterations=None,
                orbital_hvp=None, initial_energy_Eh=None, energy_Eh=None,
                scf_seconds=None, total_seconds=None, rss_mib=None)
            rows.append(row)
            if not (run / "output.out").exists():
                continue
            text = (run / "output.out").read_text()
            steps = re.search(r"VBSCF converged in\s+(\d+) iterations", text)
            trace = []
            if (run / "neo.tsv").exists():
                with (run / "neo.tsv").open() as handle:
                    trace = list(csv.DictReader(handle, delimiter="\t"))
            resource = (run / "resource.time").read_text()
            reason = re.search(r"Termination reason\s*:\s*(.+)", text)
            # The top-level exception result can print zero iterations after
            # real accepted steps. The emitted trajectory preserves their count.
            row.update(
                status="converged" if steps else "failed",
                termination_reason=reason[1] if reason else None,
                iterations=int(steps[1]) if steps else (int(trace[-1]["iteration"]) if trace else 0),
                orbital_hvp=sum(int(entry["orbital_hvp_actions"]) for entry in trace),
                initial_energy_Eh=float(re.search(r"Initial total energy\s*:\s*([-\d.]+)", text)[1]),
                energy_Eh=float(re.search(r"Final total energy\s*:\s*([-\d.]+)", text)[1]),
                scf_seconds=float(re.search(r"SCF iteration wall time\s*:\s*([\d.]+)", text)[1]),
                total_seconds=float(re.search(r"End-to-end wall time\s*:\s*([\d.]+)", text)[1]),
                rss_mib=int(re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", resource)[1])/1024,
            )
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "summarize"))
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--system", choices=("240", "LOFLEA"), default="240")
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare(args.source, args.destination, args.system)
    else:
        summarize(args.source, args.destination)
