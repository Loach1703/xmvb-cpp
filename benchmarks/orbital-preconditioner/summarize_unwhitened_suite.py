"""Compare saved Slurm runs without rerunning historical executables.

Usage: python3 summarize_unwhitened_suite.py DIAGNOSTICS_ROOT OUTPUT_CSV
Only scientifically converged references count as performance comparisons.
"""
import csv
import hashlib
from pathlib import Path
import re
import sys


def capture(pattern, text):
    match = re.search(pattern, text)
    if not match:
        raise ValueError(f'Missing field: {pattern}')
    return match.group(1)


def read_run(run):
    output = (run/'output.out').read_text()
    resources = (run/'resource.time').read_text()
    trace = list(csv.DictReader((run/'neo.tsv').open(), delimiter='\t'))
    wall = 0.0
    for part in capture(r'Elapsed.*?:\s+([\d:.]+)', resources).split(':'):
        wall = wall * 60 + float(part)
    termination = capture(r'Termination reason\s+:\s+(\S+)', output)
    exit_code = int((run/'run.exit').read_text())
    return dict(steps=len(trace), hvp=sum(int(row['orbital_hvp_actions']) for row in trace),
                energy_Eh=capture(r'Final total energy\s+:\s+(\S+)', output),
                wall_s=wall,
                application_wall_s=float(capture(r'End-to-end wall time\s+:\s+(\S+)', output)),
                rss_MiB=int(capture(r'Maximum resident set size \(kbytes\):\s+(\d+)', resources))/1024,
                node=(run/'node.txt').read_text().strip(), termination=termination,
                converged=exit_code == 0 and termination in {
                    'neo_initial_tolerance', 'neo_projected_gradient_tolerance'},
                exit_code=exit_code, run_path=str(run.resolve()))


def input_hashes(run):
    # Historical transfers sometimes omitted large basis assets. Their
    # original Slurm manifests still certify exactly which files were used.
    recorded = {}
    for line in (run/'hashes.txt').read_text().splitlines():
        digest, filename = line.split(maxsplit=1)
        name = Path(filename.strip()).name
        if name == 'input.xmi' or name.endswith('.gbs'):
            recorded[name] = digest
    for name, digest in recorded.items():
        path = run/name
        if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f'{path} no longer matches its execution manifest')
    if 'input.xmi' not in recorded:
        raise ValueError(f'{run}: no input provenance')
    return recorded


def main():
    root = Path(sys.argv[1]).resolve()
    output_csv = Path(sys.argv[2])
    rows = []
    for case in ['F2-OEO', '241', '240', '7975', 'QOCGUM', 'CERRAS', 'LOFLEA']:
        batch = '1251994' if case in {'CERRAS', 'LOFLEA'} else '1251989'
        run = root/f'raw-chart-suite-{batch}'/case
        reference_batch = ('current-suite-1249844' if case in {'F2-OEO', '240', '7975'}
                           else 'joint-isolated-1251207' if case == 'LOFLEA'
                           else 'joint-response-1251180')
        reference = root/reference_batch/case
        new = read_run(run)
        old = read_run(reference)
        same_input = input_hashes(run) == input_hashes(reference)
        comparable = same_input and new['converged'] and old['converged']
        row = dict(case=case, **new,
                   **{f'baseline_{key}': value for key, value in old.items()},
                   input_hashes_equal=same_input, performance_comparable=comparable,
                   energy_delta_Eh=float(new['energy_Eh'])-float(old['energy_Eh']),
                   wall_speedup=old['wall_s']/new['wall_s'] if comparable else '')
        rows.append(row)
    with output_csv.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print('| Case | Steps old/new | HVP old/new | Wall s old/new | RSS MiB old/new | Energy difference Eh | Paired |')
    print('|---|---:|---:|---:|---:|---:|---|')
    for row in rows:
        print(f"| {row['case']} | {row['baseline_steps']}/{row['steps']} | "
              f"{row['baseline_hvp']}/{row['hvp']} | "
              f"{row['baseline_wall_s']:.2f}/{row['wall_s']:.2f} | "
              f"{row['baseline_rss_MiB']:.1f}/{row['rss_MiB']:.1f} | "
              f"{row['energy_delta_Eh']:.2e} | {row['performance_comparable']} |")


if __name__ == '__main__':
    main()
