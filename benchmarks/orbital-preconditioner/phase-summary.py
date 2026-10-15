"""Extract paired native resources and inclusive phase totals into CSV."""

import csv
from pathlib import Path
import re
import sys


def field(pattern, text):
    found = re.search(pattern, text)
    return found.group(1) if found else ""


def read(path):
    return path.read_text() if path.exists() else ""


writer = csv.writer(sys.stdout)
writer.writerow(["case", "variant", "node", "converged", "steps", "energy_Eh",
                 "wall_s", "peak_RSS_MiB", "stage", "calls", "stage_wall_s",
                 "projector_peak_payload_bytes", "stage_cpu_s", "stage_cpu_per_wall"])
for case in sorted(Path(sys.argv[1]).iterdir()):
    if not case.is_dir():
        continue
    for variant in ("baseline", "updated"):
        run = case / variant
        output = read(run / "output.out")
        resources = read(run / "resource.time")
        diagnostics = read(run / "run.err")
        wall = field(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)", resources)
        seconds = 0.0
        for component in wall.split(":") if wall else []:
            seconds = 60 * seconds + float(component)
        rss = field(r"Maximum resident set size \(kbytes\):\s*(\d+)", resources)
        steps = field(r"VBSCF converged in\s+(\d+) iterations", output)
        common = [case.name, variant, read(case / "node.txt").strip(), bool(steps), steps,
                  field(r"Final total energy\s*:\s*([-\d.]+)", output),
                  f"{seconds:.2f}" if wall else "",
                  f"{int(rss) / 1024:.1f}" if rss else ""]
        payload = field(r"XMVB_WORKSPACE,projector_peak_payload_bytes,(\d+)", diagnostics)
        cpu = dict(re.findall(r"XMVB_CPU,(\w+),([\d.]+)", diagnostics))
        for stage, calls, duration in re.findall(r"XMVB_PHASE,(\w+),(\d+),([\d.]+)", diagnostics):
            stage_cpu = cpu.get(stage, "")
            ratio = float(stage_cpu) / float(duration) if stage_cpu and float(duration) > 0 else None
            writer.writerow([*common, stage, calls, duration, payload, stage_cpu,
                             f"{ratio:.3f}" if ratio is not None else ""])
