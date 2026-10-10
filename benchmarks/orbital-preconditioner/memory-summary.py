"""Extract exact peak heap bytes from heaptrack's Massif-compatible output."""

import csv
from pathlib import Path
import re
import sys


writer = csv.writer(sys.stdout)
writer.writerow(["case", "variant", "peak_heap_bytes", "peak_heap_MiB",
                 "peak_bytes_in_input_loading_stacks", "peak_bytes_in_chart_stacks"])
for case in sorted(Path(sys.argv[1]).iterdir()):
    if not case.is_dir():
        continue
    for variant in ("baseline", "updated"):
        heap = case / variant / "heap"
        if not (heap / "memory.massif").exists():
            continue
        peaks = [int(value) for value in re.findall(r"mem_heap_B=(\d+)",
                                                   (heap / "memory.massif").read_text())]
        peak = max(peaks) if peaks else None
        loading, chart = 0, 0
        for line in (heap / "memory-peak.stacks").read_text().splitlines():
            stack, weight = line.rsplit(" ", 1)
            if "load_vbscf_input" in stack:
                loading += int(weight)
            if "OrbitalChart" in stack:
                chart += int(weight)
        writer.writerow([case.name, variant, peak,
                         f"{peak / 1048576:.6f}" if peak is not None else "",
                         loading, chart])
