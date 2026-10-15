"""Compare saved current NEO runs with the historical CSV, without rerunning it."""

import argparse
import csv
import hashlib
import re
import statistics
from pathlib import Path

from neo_history_results import extract


def field(pattern, text):
    match = re.search(pattern, text)
    return float(match.group(1)) if match else None


def summarize(root, baseline_csv):
    with baseline_csv.open() as stream:
        baseline = {(row['case'], row['method']): row for row in csv.DictReader(stream)}
    rows = []
    for directory in sorted(root.iterdir()):
        if not directory.is_dir():
            continue
        case = directory.name
        old = baseline[(case, 'neo')]
        saved = directory / 'historical'
        with (saved / 'optimizer.tsv').open() as stream:
            old_trace = list(csv.DictReader(stream, delimiter='\t'))
        old_output = (saved / 'run.out').read_text()
        row = dict(
            case=case, orbital_type=old['orbital_type'],
            active_electrons=int(old['n_active_electrons']),
            active_orbitals=int(old['n_active_orbitals']),
            structures=int(old['n_structures']), states=int(old['n_states']),
            baseline_revision=old['revision'],
            baseline_reported_steps=int(old['accepted_steps']),
            baseline_accepted_points=sum(int(t['keyframes']) for t in old_trace),
            baseline_hvp=int(old['orbital_hvp_actions']),
            baseline_structure_actions=int(old['structure_response_actions']),
            baseline_wall_s=float(old['total_wall_s']),
            baseline_rss_mib=int(old['peak_rss_kb']) / 1024,
            baseline_energy=float(old['energy_final_hartree']),
            baseline_physical_gradient=field(
                r'Final physical \|g\|_\(G\^-1\)\s*:\s*([\d.eE+-]+)', old_output),
            baseline_input_sha256=hashlib.sha256((saved / 'input.xmi').read_bytes()).hexdigest(),
        )
        current = directory / 'current'
        status_file = current / 'status.txt'
        if not status_file.exists():
            row['status'] = 'running_or_pending'
        elif status_file.read_text().strip() != '0':
            row['status'] = 'failed'
            row['error_tail'] = (current / 'run.err').read_text()[-600:].replace('\n', ' ')
        else:
            try:
                result = extract(current)
            except ValueError as error:
                row['status'] = 'unverified'
                row['error_tail'] = str(error)
            else:
                row['status'] = 'converged'
                for key in ('reported_steps', 'accepted_points', 'hvp', 'structure_actions',
                            'wall_s', 'rss_mib', 'energy', 'physical_gradient',
                            'executable_sha256', 'input_sha256'):
                    row['current_' + key] = result[key]
                row['input_matches'] = result['input_sha256'] == row['baseline_input_sha256']
                row['energy_difference_Eh'] = result['energy'] - row['baseline_energy']
                row['wall_speedup'] = row['baseline_wall_s'] / result['wall_s']
                row['current_initial_energy'] = field(
                    r'Initial total energy\s*:\s*([-\d.]+)', (current / 'output.out').read_text())
                row['initial_energy_difference_Eh'] = (
                    row['current_initial_energy'] - float(old['energy_initial_hartree']))
                current_output = (current / 'output.out').read_text()
                row['current_projected_gradient_inf'] = field(
                    r'Final projected \|g\|_inf\s*:\s*([\d.eE+-]+)', current_output)
                row['baseline_scf_wall_s'] = float(old['scf_wall_s'])
                row['current_scf_wall_s'] = field(
                    r'SCF iteration wall time\s*:\s*([\d.eE+-]+)', current_output)
                lbfgs = baseline.get((case, 'block_lbfgs'))
                if lbfgs:
                    row['historical_lbfgs_steps'] = int(lbfgs['accepted_steps'])
                    row['historical_lbfgs_wall_s'] = float(lbfgs['total_wall_s'])
                    row['historical_lbfgs_energy'] = float(lbfgs['energy_final_hartree'])
        rows.append(row)
    return rows


def report(rows, root):
    complete = [r for r in rows if r['status'] == 'converged']
    lines = [
        '# NEO multi-system validation, 2026-10-11', '',
        'Historical reference: `results/history/paired_runs.csv`, NEO rows at '
        '`c2776c4`. Current executable: restored candidate strategy from '
        '`d86a348`, tested at `28f3a99`.', '',
        'Hanhai25 Slurm array `281497`; 22 tasks submitted concurrently. '
        'Each current run uses 32 compact-affinity OpenMP cores, one BLAS thread, '
        'exact Libcint integrals and Davidson. Thresholds: projected-gradient '
        'infinity norm 1e-3 and adjacent energy change 1e-7 Eh; existing additional '
        'OEO stationarity checks are unchanged. No optimizer code was changed.', '',
        'Inputs were copied from the saved historical run directories, not '
        'regenerated. File hashes and initial energies are compared by code. '
        'Historical timings are saved CSV measurements, not same-node replays; '
        'their affinity and system load are not controlled by this round. '
        'Timing ratios describe observed results, not an isolated causal speedup.', '',
        'Wall time is application end-to-end time. RSS is GNU-time process peak '
        'memory in MiB. Points are sums of the optimizer trace `keyframes` column '
        'for both versions: historical reported macro counts can hide multiple '
        'accepted points. Both reported counts are also retained in the CSV.', '',
        f'Converged: {len(complete)}/{len(rows)} retrieved samples. '
        f'Other statuses: {", ".join(r["case"] + ": " + r["status"] for r in rows if r["status"] != "converged") or "none"}.', '',
    ]
    if complete:
        faster = sum(r['wall_speedup'] > 1 for r in complete)
        near = sum(abs(r['energy_difference_Eh']) <= 1e-7 for r in complete)
        lines += [
            f'Current wall time is lower for {faster}/{len(complete)} completed '
            f'samples; median historical/current wall ratio is '
            f'{statistics.median(r["wall_speedup"] for r in complete):.3f}. '
            f'Final energies differ by at most 1e-7 Eh for {near}/{len(complete)}.', '',
            f'Input hashes match for {sum(r["input_matches"] for r in complete)}/'
            f'{len(complete)}; maximum absolute initial-energy difference is '
            f'{max(abs(r["initial_energy_difference_Eh"]) for r in complete):.3e} Eh.', '',
            'Passing an adjacent-step energy criterion does not guarantee final '
            'energies agree to that tolerance. In particular, historical FeCl2 '
            'has physical dual gradient 0.001149 while the current run passes '
            'its additional physical stationarity gate. The criteria are not '
            'identical across versions despite identical numerical tolerances. '
            'Physical dual norm is not the projected infinity norm; a HAO dual '
            'norm above 1e-3 alone does not violate its infinity-norm criterion.', '',
        ]
    lines += [
        '## Results', '',
        'Every paired cell is **historical / current**. Delta E is current '
        'minus historical, in Eh.', '',
        '| Case | Active e/o | States | Points | HVP | Structure actions | Wall s | RSS MiB | Delta E |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |',
    ]
    for r in rows:
        prefix = f'| {r["case"]} | {r["active_electrons"]}/{r["active_orbitals"]} | {r["states"]}'
        if r['status'] != 'converged':
            lines.append(prefix + f' | {r["status"]} | — | — | — | — | — |')
        else:
            lines.append(prefix +
                         f' | {r["baseline_accepted_points"]}/{r["current_accepted_points"]}'
                         f' | {r["baseline_hvp"]}/{r["current_hvp"]}'
                         f' | {r["baseline_structure_actions"]}/{r["current_structure_actions"]}'
                         f' | {r["baseline_wall_s"]:.3f}/{r["current_wall_s"]:.3f}'
                         f' | {r["baseline_rss_mib"]:.1f}/{r["current_rss_mib"]:.1f}'
                         f' | {r["energy_difference_Eh"]:+.3e} |')
    lines += [
        '', '## Assessment', '',
        'The MnF2 recovery is not sufficient to claim a broadly superior NEO '
        'strategy. The panel identifies both faster systems and regressions. '
        '241 needs few orbital HVPs but many more structure actions than the '
        'historical run; NEGPIA increases accepted points and orbital HVPs. '
        'These are different issues, not one molecule-specific threshold problem.', '',
        'Action counts alone do not identify the fraction of wall time spent '
        'inside a kernel. Any next optimization should measure structure '
        'certification/model refresh cost on 241 separately from the orbital '
        'trajectory regression on NEGPIA. Preserve current correctness and '
        'memory checks rather than restoring the whole historical solver.', '',
        f'Raw evidence directory: `{root}`. The comparison CSV also retains '
        'historical block-LBFGS steps, wall time and energy where available; '
        'those rows were not rerun in this round.', '',
    ]
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--baseline', type=Path, default=Path('results/history/paired_runs.csv'))
    parser.add_argument('--csv', type=Path, required=True)
    parser.add_argument('--markdown', type=Path)
    parser.add_argument('--reference', type=Path,
                        help='Previous current-version comparison CSV for matched input checks')
    args = parser.parse_args()
    rows = summarize(args.root, args.baseline)
    if args.reference:
        with args.reference.open() as stream:
            reference = {r['case']: r for r in csv.DictReader(stream)}
        for row in rows:
            if row['status'] != 'converged':
                continue
            previous = reference[row['case']]
            if previous['current_input_sha256'] != row['current_input_sha256']:
                raise ValueError(f"reference input differs for {row['case']}")
            for key in ('accepted_points', 'hvp', 'structure_actions', 'wall_s',
                        'rss_mib', 'energy', 'physical_gradient', 'executable_sha256'):
                row['reference_' + key] = previous['current_' + key]
            row['reference_wall_speedup'] = float(previous['current_wall_s']) / row['current_wall_s']
            row['energy_difference_from_reference_Eh'] = (
                row['current_energy'] - float(previous['current_energy']))
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with args.csv.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    if args.markdown:
        args.markdown.write_text(report(rows, args.root))
    print('case,status,accepted_points old/new,HVP old/new,wall old/new,RSS old/new,deltaE')
    for row in rows:
        if row['status'] != 'converged':
            print(row['case'], row['status'], row.get('error_tail', ''))
            continue
        print(f"{row['case']},{row['status']},"
              f"{row['baseline_accepted_points']}/{row['current_accepted_points']},"
              f"{row['baseline_hvp']}/{row['current_hvp']},"
              f"{row['baseline_wall_s']:.3f}/{row['current_wall_s']:.3f},"
              f"{row['baseline_rss_mib']:.1f}/{row['current_rss_mib']:.1f},"
              f"{row['energy_difference_Eh']:+.3e}")


if __name__ == '__main__':
    main()
