"""Read native run logs and print a CSV comparison; never infer convergence."""

import csv
from pathlib import Path
import re
import sys


def match(pattern, text):
    found = re.search(pattern, text)
    return found.group(1) if found else ""


writer = csv.writer(sys.stdout)
writer.writerow(["case", "iscf", "variant", "converged", "steps", "energy_Eh",
                 "wall_s", "peak_RSS_MiB", "node", "exit_code"])
for run in sorted(Path(sys.argv[1]).iterdir()):
    if not run.is_dir():
        continue
    name = re.fullmatch(r"(.+)-iscf(\d+)-(baseline|updated)", run.name)
    if not name:
        continue
    output = (run / "output.out").read_text() if (run / "output.out").exists() else ""
    resources = (run / "resource.time").read_text() if (run / "resource.time").exists() else ""
    steps = match(r"VBSCF converged in\s+(\d+) iterations", output)
    wall = match(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)", resources)
    seconds = 0.0
    for component in wall.split(":") if wall else []:
        seconds = 60 * seconds + float(component)
    rss = match(r"Maximum resident set size \(kbytes\):\s*(\d+)", resources)
    node = (run / "node.txt").read_text().strip() if (run / "node.txt").exists() else ""
    status = (run / "exit_code.txt").read_text().strip() if (run / "exit_code.txt").exists() else ""
    writer.writerow([*name.groups(), bool(steps), steps,
                     match(r"Final total energy\s*:\s*([-\d.]+)", output),
                     f"{seconds:.2f}" if wall else "",
                     f"{int(rss) / 1024:.1f}" if rss else "", node, status])
