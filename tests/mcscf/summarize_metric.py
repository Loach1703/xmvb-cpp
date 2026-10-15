#!/usr/bin/env python3
"""Parse saved metric-point audits, bounded captures and accepted-state probes."""

import argparse
import csv
from pathlib import Path
import re

from summarize_factor_gradient import save, values


def scalar(text, pattern):
    match = re.search(pattern, text)
    return float(match[1]) if match else None


def audit(artifact):
    data = values(artifact / "audit.out")
    return dict(run=artifact.name, **data,
                relative_skew=float(data["action_skew"])/float(data["action_norm"]),
                relative_feature_gap=float(data["action_feature_gap"])/float(data["feature_gram_norm"]))


def capture(artifact):
    text = (artifact / "output.out").read_text()
    with (artifact / "neo.tsv").open() as handle:
        trace = list(csv.DictReader(handle, delimiter="\t"))
    reason = re.search(r"Termination reason\s*:\s*(.+)", text)
    reason = reason[1] if reason else "unknown"
    resource = (artifact / "resource.time").read_text()
    return dict(
        run=artifact.name, reason=reason,
        application_exit=int((artifact / "run.exit").read_text()) if (artifact / "run.exit").exists()
            else int(scalar(resource,r"Exit status:\s*(\d+)")),
        accepted_steps=int(trace[-1]["iteration"]) if trace else 0,
        recorded_orbital_hvp=sum(int(row["orbital_hvp_actions"]) for row in trace),
        initial_energy_Eh=scalar(text,r"Initial total energy\s*:\s*([-\d.]+)"),
        final_energy_Eh=scalar(text,r"Final total energy\s*:\s*([-\d.]+)"),
        scf_seconds=scalar(text,r"SCF iteration wall time\s*:\s*([\d.]+)"),
        total_seconds=scalar(text,r"End-to-end wall time\s*:\s*([\d.]+)"),
        rss_mib=scalar(resource,r"Maximum resident set size \(kbytes\):\s*(\d+)")/1024,
        raw_orbitals_saved=(artifact / "orbitals.bin").is_file())


def accepted(artifact):
    return dict(run=artifact.name, kind="synthetic_production_gradient",
                **values(artifact / "accepted-state16.out"))


def hvp(artifact):
    return dict(run=artifact.name, kind="synthetic_production_hvp",
                **values(artifact / "accepted-hvp16.out"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("audit", "capture", "accepted", "hvp"))
    parser.add_argument("output", type=Path)
    parser.add_argument("artifacts", type=Path, nargs="+")
    args = parser.parse_args()
    extract = dict(audit=audit, capture=capture, accepted=accepted, hvp=hvp)[args.kind]
    save(args.output,[extract(artifact) for artifact in args.artifacts])
