"""Summarize scalar NEO diagnostics; no manual transcription of trajectories."""

import argparse
import collections
import csv
import json
from pathlib import Path


def summarize(path):
    trials = collections.defaultdict(lambda: collections.defaultdict(list))
    with path.open() as stream:
        for line in stream:
            if not line.startswith("NEO_DIAG,"):
                continue
            _, keyframe, attempt, event, *fields = next(csv.reader([line]))
            values = dict(field.split("=", 1) for field in fields)
            trials[int(keyframe), int(attempt)][event].append(values)
    rows = []
    for (keyframe, attempt), events in sorted(trials.items()):
        trial = events["trial"]
        if not trial:
            continue
        row = {"keyframe": keyframe, "attempt": attempt, **trial[-1]}
        candidates = events["candidate"]
        if candidates:
            row.update({f"final_{key}": candidates[-1][key]
                        for key in ("ro", "rs", "target", "no", "ns")})
        basis = events["basis"]
        independent = [x for x in basis if x["independent"] == "1"]
        row["orbital_sources"] = dict(collections.Counter(
            x["source"] for x in independent))
        row["dependent_proposals"] = sum(x["independent"] == "0" for x in basis)
        row["min_retained_fraction"] = min(
            (float(x["retained"]) for x in independent), default=None)
        pcg = events["precondition"]
        row["pcg_iterations"] = sum(int(x["iterations"]) for x in pcg)
        row["pcg_calls"] = len(pcg)
        checks = events["response_check"]
        row["response_checks"] = len(checks)
        row["response_reexpansions"] = sum(x["retained"] == "1" for x in checks)
        row["response_max_norm_ratio"] = max(
            (float(x["after"]) / float(x["before"])
             for x in checks if float(x["before"]) > 0), default=None)
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args()
    rows = summarize(args.log)
    if args.csv and rows:
        with args.csv.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            for row in rows:
                writer.writerow({**row, "orbital_sources": json.dumps(
                    row["orbital_sources"], sort_keys=True)})
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
