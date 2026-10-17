"""Collect exact-support growth and explicitly sampled connection estimates."""
import csv
import json
from pathlib import Path
import re

root=Path(__file__).resolve().parent
rows=[]
for p in sorted((root/'completed/run-1257163').glob('*/result.json')):
    d=json.loads(p.read_text())
    c=json.loads((p.parent/'connectivity.json').read_text())
    text=(p.parent/'resource.time').read_text()
    rss=int(re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',text)[1])
    rows.append(dict(n=d['n'],initial_support=2**(d['n']//2),
                     transformed_support=d['support'],transform_seconds=d['transform_seconds'],
                     coefficients_mib=d['coefficient_bytes']/2**20,peak_rss_mib=rss/1024,
                     connection_sample_rows=c['sampled_rows'],
                     mean_potential_connections=c['mean_potential_connections_in_support'],
                     mean_standard_error=c['mean_standard_error'],
                     estimated_directed_connections=c['estimated_directed_connections']))
with (root/'summary.csv').open('w',newline='') as handle:
    writer=csv.DictWriter(handle,fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
