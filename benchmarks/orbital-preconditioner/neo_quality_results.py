"""Compare diagnostic nonlinear candidates with the final trial at the same point/radius."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_events(path):
    events = defaultdict(lambda: defaultdict(list))
    with path.open() as stream:
        for line in stream:
            if not line.startswith("NEO_DIAG,"):
                continue
            _, keyframe, attempt, event, *fields = next(csv.reader([line]))
            events[int(keyframe), int(attempt)][event].append(
                dict(field.split("=", 1) for field in fields))
    return events


def candidate_rows(path, trace=None):
    gradients = {}
    if trace is not None:
        with trace.open() as stream:
            gradients = {int(row["iteration"]): float(row["accepted_physical_gradient_norm"])
                         for row in csv.DictReader(stream, delimiter="\t")}
    for (keyframe, attempt), events in sorted(read_events(path).items()):
        if not events["trial"]:
            continue
        final = events["trial"][-1]
        for sample in events["quality"]:
            row = {"keyframe": keyframe, "attempt": attempt, **sample}
            row.update({f"final_{key}": final[key]
                        for key in ("hvp", "actual", "predicted", "rho")})
            row["final_energy_acceptable"] = final["accepted"]
            accepted = final["accepted"] == "1" and not (
                events["keyframe"] and events["keyframe"][-1]["reliable"] == "0")
            row["final_accepted"] = int(accepted)
            row["final_gradient"] = gradients.get(keyframe, "") if accepted else ""
            actual = float(final["actual"])
            row["actual_fraction_of_final"] = (
                float(sample["actual"]) / actual if actual > 0 else "")
            source = float(sample["source_gradient"])
            row["gradient_ratio"] = (
                float(sample["actual_gradient"]) / source if source > 0 else "")
            yield row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--trace", type=Path, help="Macro TSV with accepted physical gradients")
    args = parser.parse_args()
    rows = list(candidate_rows(args.log, args.trace))
    if not rows:
        raise SystemExit("No nonlinear candidate audit rows found")
    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {len(rows)} observations to {args.csv}")
    print(f"Nonlinear audit time: {sum(float(x['audit_seconds']) for x in rows):.3f} s")
    for row in rows:
        if int(row["final_hvp"]) < 50:
            continue
        fraction = row["actual_fraction_of_final"]
        fraction_text = f"{fraction:.3f}" if fraction != "" else "n/a"
        print(f"keyframe={row['keyframe']} attempt={row['attempt']} "
              f"HVP={row['hvp']}/{row['final_hvp']} "
              f"energy_acceptable={row['energy_acceptable']} "
              f"actual/final={fraction_text} "
              f"gradient/source={row['gradient_ratio']:.3f} "
              f"residual={float(row['residual']):.3e}/"
              f"{float(row['target']):.3e}")


if __name__ == "__main__":
    main()
