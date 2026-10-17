#!/usr/bin/env python3
"""Prepare fixed C12 structure spaces and parse native convergence records."""

import argparse
import csv
import hashlib
import json
from math import comb
from pathlib import Path
import random
import re
import shutil
import sys

from summarize_scf_step import elapsed_seconds, field

SEED = 20261017
LABELS = ("ion01", "ion02", "full", "mixed3", "mixed20", "mixed600")


def pairings(sites):
    """Noncrossing singlet pairings form an independent Rumer basis."""
    if not sites:
        return [()]
    result = []
    for partner in range(1, len(sites), 2):
        for inside in pairings(sites[1:partner]):
            for outside in pairings(sites[partner + 1:]):
                result.append(((sites[0], sites[partner]),) + inside + outside)
    return result


def mixed_structures(count=20):
    """Select before calculation; retain the same nested mixed subsets."""
    rng = random.Random(SEED)
    sites = tuple(range(29, 41))
    selected = [tuple(zip(sites[::2], sites[1::2]))]
    while len(selected) < count:
        # The first three contain orders 0,1,2; the full set also contains 3.
        # Preserve the archived first twenty columns. Larger arbitrary sets
        # sample occupations independently, avoiding exhaustion of the 132
        # distinct covalent columns in a forced round-robin schedule.
        order = len(selected) % 4 if len(selected) < 20 else rng.randrange(4)
        occupied_and_empty = rng.sample(sites, 2 * order)
        doubly = occupied_and_empty[:order]
        singles = tuple(p for p in sites if p not in occupied_and_empty)
        bonds = tuple(sorted(tuple((p, p) for p in doubly) + rng.choice(pairings(singles))))
        if bonds in selected:
            continue
        if len(bonds) != 6:
            raise ValueError("mixed structure does not contain twelve active electrons")
        selected.append(bonds)
    return ["1:28 " + " ".join(str(p) for bond in bonds for p in bond) for bonds in selected]


def sector_size(max_order):
    return sum(comb(12, d) * comb(12 - d, d) *
               comb(12 - 2 * d, 6 - d) // (7 - d) for d in range(max_order + 1))


def prepare(source, bundle):
    original = source.read_text()
    if re.search(r"^\s*\$STR\s*$", original, re.I | re.M):
        raise ValueError("source already contains an explicit structure block")
    text, changes = re.subn(r"^iscf=.*$", "iscf=9", original, flags=re.I | re.M)
    if changes != 1:
        raise ValueError("expected one ISCF control")
    explicit = re.sub(r"^str=.*\n", "", text, flags=re.I | re.M)
    structures = mixed_structures(600)
    manifest = dict(source_input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    selection_seed=SEED, cases=[])
    for label, keyword, count in (("ion01", "ION(0:1)", sector_size(1)),
                                  ("ion02", "ION(0:2)", sector_size(2)),
                                  ("full", "FULL", sector_size(6)),
                                  ("mixed3", None, 3), ("mixed20", None, 20),
                                  ("mixed600", None, 600)):
        if keyword:
            deck, changes = re.subn(r"^str=.*$", f"str={keyword}", text, flags=re.I | re.M)
            if changes != 1:
                raise ValueError("expected one generated-structure control")
        else:
            deck = explicit + "\n$STR\n" + "\n".join(structures[:count]) + "\n$END\n"
        run = bundle / "cases" / label
        run.mkdir(parents=True, exist_ok=True)
        (run / "input.xmi").write_text(deck)
        shutil.copyfile(bundle / "cc-pvdz.gbs", run / "cc-pvdz.gbs")
        manifest["cases"].append(dict(label=label, keyword=keyword, expected_structures=count,
            input_sha256=hashlib.sha256(deck.encode()).hexdigest(),
            explicit_structures=structures[:count] if keyword is None else None))
    (bundle / "cases.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


def summarize(root):
    # FULL can print hundreds of MiB of structure weights. Read only scalar
    # diagnostics, never retain those coefficient tables in the report parser.
    labels = ("Number of structures:", "Termination reason", "Initial total energy",
              "Final total energy", "End-to-end wall time", "SCF iteration wall time")
    for label in LABELS:
        run = root / label
        text = ""
        if (run / "output.out").exists():
            with (run / "output.out").open() as stream:
                text = "".join(line.lstrip() for line in stream if line.lstrip().startswith(labels))
        resource = (run / "resource.time").read_text() if (run / "resource.time").exists() else ""
        trace = []
        if (run / "ciah.tsv").exists():
            with (run / "ciah.tsv").open() as stream:
                trace = list(csv.DictReader(stream, delimiter="\t"))
        reason = field(text, r"^Termination reason\s*:\s*(\S+)")
        done = (run / "run.exit").exists()
        exit_code = (run / "run.exit").read_text().strip() if done else None
        rss = field(resource, r"Maximum resident set size \(kbytes\):\s*(\d+)", int)
        accepted_hvp = sum(int(r["orbital_hvp_actions"]) for r in trace)
        errors = (run / "run.err").read_text() if (run / "run.err").exists() else ""
        failed_hvp = sum(map(int, re.findall(r"^CIAH_FAILURE,[^\n]*,hvp=(\d+)", errors, re.M)))
        yield dict(space=label,
            status="converged" if done and exit_code == "0" and reason == "ciah_projected_gradient_tolerance" else
                   "finished_without_convergence" if done else "running" if text else "pending",
            termination=reason, exit_code=exit_code,
            structures=field(text, r"Number of structures:\s*(\d+)", int),
            steps=len(trace), hvp=accepted_hvp + failed_hvp,
            accepted_hvp=accepted_hvp, failed_hvp=failed_hvp,
            responses=sum(int(r["structure_response_actions"]) for r in trace),
            initial_Eh=field(text, r"^Initial total energy\s*:\s*([-\d.]+)"),
            final_Eh=field(text, r"^Final total energy\s*:\s*([-\d.]+)"),
            application_s=field(text, r"^End-to-end wall time\s*:\s*([\d.]+)", float),
            scf_s=field(text, r"^SCF iteration wall time\s*:\s*([\d.]+)", float),
            process_s=elapsed_seconds(resource), rss_GiB=rss / 1048576 if rss is not None else None,
            node=(run / "node.txt").read_text().strip() if (run / "node.txt").exists() else None)


def compare(control, current):
    """Compare identical inputs; partial records never produce speedups."""
    baseline = {row["space"]: row for row in summarize(control)}
    for row in summarize(current):
        old = baseline[row["space"]]
        paths = [root / row["space"] / "input.xmi" for root in (control, current)]
        if all(path.exists() for path in paths):
            if paths[0].read_bytes() != paths[1].read_bytes():
                raise ValueError(f"different inputs in {row['space']}")
        matched = old["status"] == row["status"] == "converged"
        result = dict(space=row["space"], structures=row["structures"])
        for name in ("status", "steps", "hvp", "final_Eh", "application_s", "process_s", "rss_GiB"):
            result[f"control_{name}"] = old[name]
            result[f"current_{name}"] = row[name]
        result["energy_delta_Eh"] = float(row["final_Eh"]) - float(old["final_Eh"]) if matched else None
        result["wall_speedup"] = old["application_s"] / row["application_s"] if matched else None
        yield result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("source", type=Path)
    prep.add_argument("bundle", type=Path)
    summary = sub.add_parser("summarize")
    summary.add_argument("root", type=Path)
    comparison = sub.add_parser("compare")
    comparison.add_argument("control", type=Path)
    comparison.add_argument("current", type=Path)
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare(args.source, args.bundle)
    else:
        rows = list(compare(args.control, args.current) if args.mode == "compare" else summarize(args.root))
        writer = csv.DictWriter(sys.stdout, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
