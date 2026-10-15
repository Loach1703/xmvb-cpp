"""Compare new NEO frame-only runs with saved unwhitened regressions."""

import argparse
import csv
from pathlib import Path

from summarize_unwhitened_suite import input_hashes, read_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with args.baseline.open() as stream:
        references = {row["case"]: row for row in csv.DictReader(stream)}
    rows = []
    for case in ("FeCl2", "MnF2", "241", "240", "F2-OEO", "F2-SA2"):
        old = references[case]
        run = args.runs / case
        new = read_run(run)
        same_input = input_hashes(run) == input_hashes(Path(old["run_path"]))
        row = dict(case=case, **new, input_hashes_equal=same_input,
                   baseline_steps=int(old["steps"]), baseline_hvp=int(old["hvp"]),
                   baseline_energy_Eh=old["energy_Eh"], baseline_wall_s=float(old["wall_s"]),
                   baseline_rss_MiB=float(old["rss_MiB"]), baseline_node=old["node"],
                   energy_delta_Eh=float(new["energy_Eh"])-float(old["energy_Eh"]),
                   binary_sha256=(run/"hashes.txt").read_text().split()[0])
        if not same_input or not new["converged"]:
            raise ValueError(f"{case}: not a matched converged comparison")
        rows.append(row)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print("| Case | Steps old/new | HVP old/new | Wall s old/new | RSS MiB old/new | Delta E / Eh |")
    print("|---|---:|---:|---:|---:|---:|")
    for row in rows:
        print(f"| {row['case']} | {row['baseline_steps']}/{row['steps']} | "
              f"{row['baseline_hvp']}/{row['hvp']} | "
              f"{row['baseline_wall_s']:.2f}/{row['wall_s']:.2f} | "
              f"{row['baseline_rss_MiB']:.1f}/{row['rss_MiB']:.1f} | "
              f"{row['energy_delta_Eh']:.2e} |")


if __name__ == "__main__":
    main()
