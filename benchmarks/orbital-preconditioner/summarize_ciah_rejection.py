"""Extract bounded late-point CIAH diagnostics; do not label them converged."""

import argparse
import csv
from pathlib import Path
import re

from summarize_unwhitened_suite import read_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = []
    for spin in ("singlet", "triplet"):
        run = args.runs / spin
        metrics = read_run(run)
        events = (run / "run.err").read_text()
        trials = re.findall(r"ciah_macro_trial,[^\n]+", events)
        accepted = [event for event in trials if event.endswith("accepted=1")]
        if len(accepted) != metrics["steps"]:
            raise ValueError(f"{spin}: accepted trace and event counts disagree")
        origin = float(re.search(r"current_energy=([^,]+)", trials[0])[1])
        ci_trials = [float(value) for value in re.findall(
            r"XMVB_PHASE_EVENT,ciah_trial,end,([^\n]+)", events)]
        row = dict(spin=spin, seed="saved-late-point-molden", **metrics,
                   initial_energy_Eh=origin,
                   energy_decrease_mEh=1000 * (origin-float(metrics["energy_Eh"])),
                   energy_model_stops=events.count("ciah_keyframe_energy_model,"),
                   origin_restarts=events.count("ciah_origin_restart,"),
                   trial_count=len(ci_trials), trial_wall_s=sum(ci_trials),
                   binary_sha256=(run / "hashes.txt").read_text().split()[0])
        rows.append(row)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print("| Spin | Steps | HVP | Wall s | RSS GiB | Energy / Eh | Status |")
    print("|---|---:|---:|---:|---:|---:|---|")
    for row in rows:
        print(f"| {row['spin']} | {row['steps']} | {row['hvp']} | "
              f"{row['wall_s']:.2f} | {row['rss_MiB']/1024:.2f} | "
              f"{float(row['energy_Eh']):.12f} | {row['termination']} |")


if __name__ == "__main__":
    main()
