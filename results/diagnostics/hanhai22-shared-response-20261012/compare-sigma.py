"""Compare saved runs with the exterior-transform baseline; no workloads run.

Usage: python3 compare-sigma.py sigma-batched-1250771 sigma-reduction-1250785
Output: CSV on stdout. Raw artifacts are kept beside this script.
"""
import csv
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parent
BASELINE = ROOT / "transform-clean-1250728"


def read_run(run):
    output = (run / "output.out").read_text()
    match = re.search(r"VBSCF converged in\s+(\d+) iterations", output)
    if (run / "run.exit").read_text().strip() != "0" or match is None:
        raise ValueError(f"run did not converge: {run}")
    with (run / "neo.tsv").open() as trace:
        hvp = sum(int(row["orbital_hvp_actions"])
                  for row in csv.DictReader(trace, delimiter="\t"))
    memory = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)",
                       (run / "resource.time").read_text())
    inputs = {Path(parts[-1]).name: parts[0]
              for parts in map(str.split, (run / "hashes.txt").read_text().splitlines())
              if parts and (parts[-1].endswith(".gbs") or parts[-1].endswith(".xmi"))}
    return dict(steps=int(match[1]), hvp=hvp,
                energy_Eh=float(re.search(r"Final total energy\s*:\s*([-\d.]+)", output)[1]),
                wall_s=float(re.search(r"End-to-end wall time\s*:\s*([\d.]+)", output)[1]),
                rss_MiB=int(memory[1]) / 1024,
                node=(run / "node.txt").read_text().strip(), inputs=inputs)


writer = csv.writer(sys.stdout, lineterminator="\n")
fields = ("steps", "hvp", "energy_Eh", "wall_s", "rss_MiB", "node")
writer.writerow(("run", "case", *fields, *("baseline_" + f for f in fields),
                 "energy_delta_Eh", "speedup", "input_hashes_equal"))
for name in sys.argv[1:]:
    for exit_file in sorted((ROOT / name).glob("*/run.exit")):
        run = exit_file.parent
        before = read_run(BASELINE / run.name)
        after = read_run(run)
        writer.writerow((name, run.name, *(after[f] for f in fields),
                         *(before[f] for f in fields),
                         after["energy_Eh"] - before["energy_Eh"],
                         before["wall_s"] / after["wall_s"],
                         after["inputs"] == before["inputs"]))
