"""Extract completed production-default regressions and their native resources."""
import csv
from pathlib import Path
import re

root = Path(__file__).resolve().parent
destination = root.parents[2]/'benchmarks/orbital-preconditioner/neo-unwhitened-default-20261013.csv'
rows = []
for batch in ['1252070', '1252077', '1252083']:
    for run in sorted((root/f'neo-raw-regression-{batch}').iterdir()):
        if not run.is_dir():
            continue
        report = (run/'output.out').read_text()
        time = (run/'resource.time').read_text()
        def extract(pattern, text=report):
            match = re.search(pattern, text)
            if not match:
                raise ValueError(f'{run}: missing {pattern}')
            return match.group(1)
        termination = extract(r'Termination reason\s+:\s+(\S+)')
        if int((run/'run.exit').read_text()) != 0 or 'tolerance' not in termination:
            raise ValueError(f'{run}: calculation not converged')
        backend = 'tnhvp' if run.name == 'F2-TN' else 'block_lbfgs' if run.name == 'F2-block' else 'neo'
        steps = int(extract(r'VBSCF converged in\s+(\d+) iterations'))
        hvp = ''  # TN prints block counts, not the NEO direction-count convention.
        if backend == 'neo':
            if extract(r'Orbital coordinates\s+:\s+([^\n]+)').strip() != 'unwhitened quotient':
                raise ValueError(f'{run}: NEO did not select unwhitened coordinates')
            trace = list(csv.DictReader((run/'neo.tsv').open(), delimiter='\t'))
            if len(trace) != steps:
                raise ValueError(f'{run}: inconsistent accepted iteration counts')
            hvp = sum(int(row['orbital_hvp_actions']) for row in trace)
        elif backend == 'block_lbfgs':
            hvp = 0
        wall = 0.0
        for part in extract(r'Elapsed.*?:\s+([\d:.]+)', time).split(':'):
            wall = wall*60 + float(part)
        rows.append(dict(case=run.name, backend=backend, steps=steps, hvp=hvp,
                         energy_Eh=extract(r'Final total energy\s+:\s+(\S+)'),
                         wall_s=wall,
                         rss_MiB=int(extract(r'Maximum resident set size \(kbytes\):\s+(\d+)', time))/1024,
                         threads=32, node=(run/'node.txt').read_text().strip(),
                         binary_sha256=(run/'hashes.txt').read_text().split()[0],
                         termination=termination, run_path=str(run)))
with destination.open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
for row in rows:
    print(row)
