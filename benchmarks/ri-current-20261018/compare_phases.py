#!/usr/bin/env python3
"""Compare saved phase timers (nested scopes, not additive totals)."""

import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate_log", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    timers = {}
    for line in args.candidate_log.read_text().splitlines():
        if line.startswith("XMVB_PHASE,"):
            _, phase, calls, seconds = line.split(",")
            if calls == "calls":
                continue
            timers[phase] = (int(calls), float(seconds))
    with args.baseline.open() as stream:
        rows = list(csv.DictReader(stream))
    recorded = {row["phase"] for row in rows}
    rows.extend(dict(phase=phase,exact_calls="",exact_s="",ri_calls="",ri_s="")
                for phase in timers if phase not in recorded)
    for row in rows:
        calls, seconds = timers[row["phase"]]
        row["candidate_calls"], row["candidate_s"] = calls, seconds
        row["same_call_count"] = int(row["ri_calls"]) == calls if row["ri_calls"] else ""
        row["ri_over_candidate_time"] = float(row["ri_s"])/seconds if seconds and row["ri_s"] else ""
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
