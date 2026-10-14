"""Compare saved singlet traces; no electronic-structure calculations."""
import csv
import gzip
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent


def read_run(name):
    folder = ROOT / name / "singlet"
    output_path = folder / "output.out"
    output = output_path.read_text() if output_path.exists() else gzip.open(
        str(output_path)+".gz", "rt").read()
    energies = {int(i): float(e) for i, e in re.findall(
        r"^\s+(\d+)\s+(-2244\.\d+)\s+[-\d.]+\s+[\d.Ee+-]+", output, re.M)}
    path = folder / "run.err"
    stream = path.open() if path.exists() else gzip.open(str(path)+".gz", "rt")
    events = {}
    with stream:
        for line in stream:
            f = line.strip().split(",")
            if len(f) < 5 or f[0] != "NEO_DIAG":
                continue
            events.setdefault(int(f[1]), []).append((f[3], dict(
                x.split("=", 1) for x in f[4:] if "=" in x)))
    rows = []
    cumulative = 0
    for row in csv.DictReader((folder / "neo.tsv").open(), delimiter="\t"):
        i = int(row["iteration"])
        hvp = int(row["orbital_hvp_actions"])
        cumulative += hvp
        selected = events.get(i, [])
        stop = next((v["reason"] for k, v in selected if k == "ciah_micro_end"), "")
        interpolation = [v for k, v in selected if k == "ciah_segment_interpolation"]
        rows.append(dict(run=name, macro=i, energy=energies[i], hvp=hvp,
                         cumulative_hvp=cumulative, keyframes=int(row["keyframes"]),
                         actual_reduction=float(row["actual_reduction"]),
                         physical_gradient=float(row["accepted_physical_gradient_norm"]),
                         stop=stop, interpolations=interpolation))
    return rows


def main():
    before = read_run("porphine-horizontal-1254797")
    after = read_run("porphine-segment-1255029")
    summary = dict(before=before, after=after,
                   hvp_through_17=[sum(r["hvp"] for r in run[:17]) for run in (before, after)],
                   additional_tail_hvp=sum(r["hvp"] for r in after[17:]),
                   first_new_macro_below_old_terminal=next(
                       r["macro"] for r in after if r["energy"] <= before[-1]["energy"]))
    target = ROOT / "porphine-segment-1255029"
    (target / "trajectory-comparison.json").write_text(json.dumps(summary, indent=2)+"\n")
    with (target / "trajectory-comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=before[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(dict(r, interpolations=json.dumps(r["interpolations"]))
                         for r in before+after)
    print(json.dumps({k: v for k, v in summary.items() if k not in ("before", "after")}))


if __name__ == "__main__":
    main()
