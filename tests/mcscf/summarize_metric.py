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
    text = (artifact / "accepted-hvp16.err").read_text()
    entries = re.findall(r"^XMVB_FACTOR_ADJOINT,([^\n]+)", text, re.M)
    # This fixture evaluates one accepted gradient, one HVP and two difference
    # gradients. Only the second integral adjoint belongs to the HVP.
    if entries and len(entries) != 4:
        raise ValueError(f"{artifact}: expected four integral-adjoint records")
    adjoint = dict(item.split("=", 1) for item in entries[1].split(",")) if entries else {}
    fields = ("bra", "ket", "left_tile", "right_tile", "channel_tile",
              "workspace_bytes", "source_loads", "load", "scatter", "covectors", "gram")
    data = values(artifact / "accepted-hvp16.out")
    reuses = data.pop("normalization_adjoint_reuses", None)
    residuals = re.findall(r"^XMVB_FACTOR_GENERATOR,hamiltonian=1,", text, re.M)
    if entries and len(residuals) < 3:
        raise ValueError(f"{artifact}: missing accepted/difference residual sweeps")
    return dict(run=artifact.name, kind="synthetic_production_hvp",
                **data, normalization_adjoint_reuses=reuses,
                hvp_residual_generator_sweeps=len(residuals)-3 if entries else None,
                **{f"integral_{key}": adjoint.get(key) for key in fields})


def tiles(artifact):
    """Generated-reader kernel proof, separate from converged SCF timing."""
    native = values(artifact / "tile-adjoint12.out")
    oracle = values(artifact / "tile-adjoint12-reference.out")
    records = re.findall(r"^XMVB_TILE,([^\n]+)",
                         (artifact / "tile-adjoint12.err").read_text(), re.M)
    if len(records) != 1:
        raise ValueError(f"{artifact}: expected one native generated channel sweep")
    work = dict(item.split("=",1) for item in records[0].split(","))
    return dict(run=artifact.name, kind="generated_support_sigma_adjoint",
                orbitals=int(native["orbitals"]), alpha_strings=int(native["alpha_strings"]),
                beta_strings=int(native["beta_strings"]),
                generated_targets=int(work["targets"]),
                full_targets=int(native["alpha_strings"])*int(native["beta_strings"]),
                batches=int(work["batches"]), peak_payload_bytes=int(work["peak_payload"]),
                generated_elements=int(work["source_elements"]),
                prepare_seconds=float(work["prepare"]), read_seconds=float(work["read"]),
                compute_seconds=float(work["compute"]), consume_seconds=float(work["write"]),
                native_seconds=float(native["sigma_adjoint_seconds"]),
                native_rss_mib=float(native["max_rss_kib"])/1024,
                oracle_dense_seconds=float(oracle["reference_seconds"]),
                oracle_process_rss_mib=float(oracle["max_rss_kib"])/1024,
                sigma_error=float(oracle["sigma_error"]),
                one_electron_error=float(oracle["one_electron_error"]),
                pair_kernel_error=float(oracle["pair_kernel_error"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("audit", "capture", "accepted", "hvp", "tiles"))
    parser.add_argument("output", type=Path)
    parser.add_argument("artifacts", type=Path, nargs="+")
    args = parser.parse_args()
    extract = dict(audit=audit, capture=capture, accepted=accepted, hvp=hvp, tiles=tiles)[args.kind]
    save(args.output,[extract(artifact) for artifact in args.artifacts])
