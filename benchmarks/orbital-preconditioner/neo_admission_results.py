#!/usr/bin/env python3
"""Write native admission-ablation results to CSV, retaining incomplete runs."""

import argparse
import csv
import re
import sys
from pathlib import Path


def extract(text, pattern):
    match = re.search(pattern, text)
    return match.group(1) if match else ""


def read(path):
    return path.read_text() if path.exists() else ""


def summarize(case, directory):
    output = read(directory / "output.out")
    resources = read(directory / "resource.time")
    trace = directory / "trace.tsv"
    if trace.exists():
        with trace.open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
    else:
        rows = []
    status = read(directory / "status.txt").strip()
    if not status:
        status = extract(read(directory / "run.err"), r"Exited with exit code (\d+)")
    converged = status == "0" and "VBSCF converged" in output
    hashes = read(directory / "hashes.txt").splitlines()
    rss = extract(resources, r"Maximum resident set size \(kbytes\): (\d+)")
    last = rows[-1] if rows else {}
    result = {
        "case": case,
        "variant": directory.name,
        "node": read(directory.parent / "node.txt").strip(),
        "status": "converged" if converged else "finished_unconverged" if status else "incomplete",
        "exit_code": status,
        "points": sum(int(row["keyframes"]) for row in rows),
        "hvp": sum(int(row["orbital_hvp_actions"]) for row in rows),
        "structure_actions": sum(int(row["structure_response_actions"]) for row in rows),
        "rejected": sum(int(row["rejected_trials"]) for row in rows),
        "total_s": extract(output, r"End-to-end wall time\s*:\s*([\d.]+)"),
        "scf_s": extract(output, r"SCF iteration wall time\s*:\s*([\d.]+)"),
        "rss_mib": float(rss) / 1024 if rss else "",
        "initial_energy": extract(output, r"Initial total energy\s*:\s*([-\d.]+)"),
        "final_energy": extract(output, r"Final total energy\s*:\s*([-\d.]+)"),
        "final_g_inf": extract(output, r"Final gradient \|g\|_inf\s*:\s*([\deE+.-]+)"),
        "last_dE": last.get("actual_reduction", ""),
        "last_physical_gradient": last.get("accepted_physical_gradient_norm", ""),
        "last_radius": last.get("accepted_trial_radius", ""),
        "input_sha256": hashes[1].split()[0] if len(hashes) > 1 else "",
        "executable_sha256": hashes[0].split()[0] if hashes else "",
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("roots", nargs="+", type=Path)
    args = parser.parse_args()
    results = []
    for root in args.roots:
        cases = [("csd_NEGPIA_vb", root)] if (root / "current").is_dir() else [
            (directory.name, directory) for directory in sorted(root.iterdir()) if directory.is_dir()
        ]
        for case, directory in cases:
            for variant in sorted(directory.iterdir()):
                if variant.is_dir() and (variant / "input.xmi").exists():
                    results.append(summarize(case, variant))
    if not results:
        parser.error("no run directories found")
    writer = csv.DictWriter(sys.stdout, fieldnames=list(results[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(results)


if __name__ == "__main__":
    main()
