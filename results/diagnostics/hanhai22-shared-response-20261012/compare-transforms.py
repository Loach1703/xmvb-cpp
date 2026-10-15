"""Parse saved native runs; never execute a chemistry workload.

Usage: python3 compare-transforms.py > transform-comparison.csv
"""
import csv
from pathlib import Path
import re
import sys

root = Path(__file__).resolve().parent
baseline = {row["case"]: row for row in csv.DictReader(
    (root / "keyframe-lean-comparison.csv").open())}
groups = {
    "transform_only": ["transform-only-1250698", "transform-only-1250709"],
    "transform_and_reuse": ["transform-clean-1250728"],
}
writer = csv.writer(sys.stdout, lineterminator="\n")
writer.writerow(("variant", "case", "converged", "steps", "hvp", "energy_Eh",
                 "wall_s", "rss_MiB", "node", "baseline_wall_s", "baseline_hvp",
                 "baseline_rss_MiB", "energy_delta_Eh", "speedup", "input_hash_equal"))

def input_hash(path):
    for line in path.read_text().splitlines():
        fields = line.split()
        if Path(fields[-1]).name == "input.xmi":
            return fields[0]
    raise ValueError(f"missing input hash in {path}")

for variant, directories in groups.items():
    for directory in directories:
        for trace in sorted((root / directory).glob("*/neo.tsv")):
            run = trace.parent
            case = run.name
            if case not in baseline or not (run / "run.exit").exists():
                continue
            previous = baseline[case]
            output = (run / "output.out").read_text()
            energy = float(re.search(r"Final total energy\s*:\s*([-\d.]+)", output)[1])
            wall = float(re.search(r"End-to-end wall time\s*:\s*([\d.]+)", output)[1])
            memory = float(re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)",
                                     (run / "resource.time").read_text())[1]) / 1024
            rows = list(csv.DictReader(trace.open(), delimiter="\t"))
            hvp = sum(int(row["orbital_hvp_actions"]) for row in rows)
            # Application exit plus explicit final convergence message.
            converged = ((run / "run.exit").read_text().strip() == "0" and
                         re.search(r"VBSCF converged in\s+\d+ iterations", output) is not None)
            same_input = input_hash(run / "hashes.txt") == input_hash(
                root / "keyframe-lean-1250648" / case / "hashes.txt")
            writer.writerow((variant, case, converged, len(rows), hvp, f"{energy:.12f}",
                             wall, memory, (run / "node.txt").read_text().strip(),
                             previous["wall_s"], previous["hvp"], previous["rss_MiB"],
                             energy - float(previous["energy"]),
                             float(previous["wall_s"]) / wall, same_input))
