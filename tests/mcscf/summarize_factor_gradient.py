#!/usr/bin/env python3
"""Parse factor-adjoint kernels and converged 240 runs; never execute workloads."""

import argparse
import csv
from itertools import product
from pathlib import Path
import re


def values(path):
    return dict(line.split("=", 1) for line in path.read_text().splitlines() if "=" in line)


def save(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def kernels(artifacts, output):
    rows = []
    for artifact in artifacts:
        for consumer, n in product(("factor-adjoint", "factor-direction"), (12, 14)):
            path = artifact / f"{n}-{consumer}.out"
            if not path.exists():
                continue
            native = values(path)
            reference = values(artifact / f"{n}-{consumer}-reference.out")
            rows.append(dict(
                run=artifact.name, consumer=consumer, orbitals=n, node=(artifact / "node.txt").read_text().strip(),
                bra_sources=int(native["bra_sources"]), ket_sources=int(native["ket_sources"]),
                cold_seconds=float(native["cold_adjoint_seconds"]),
                warm_seconds=float(native["warm_adjoint_seconds"]),
                rss_mib=float(native["native_max_rss_kib"]) / 1024,
                reference_process_factor_seconds=float(reference["warm_adjoint_seconds"]),
                product_reference_seconds=float(reference["reference_seconds"]),
                one_electron_error=float(reference["one_electron_error"]),
                pair_kernel_error=float(reference["pair_kernel_error"])))
    save(output, rows)


def scf(artifact, baseline, output):
    with baseline.open() as handle:
        old = {row["space"]: row for row in csv.DictReader(handle)
               if row["run"] == "240-metric-1260422"}
    rows = []
    for label in ("cov", "ion01", "ion02", "ion03", "full"):
        run = artifact / label
        text = (run / "output.out").read_text()
        iterations = re.search(r"VBSCF converged in\s+(\d+) iterations", text)
        if iterations is None:
            raise ValueError(f"{run}: no convergence record")
        with (run / "neo.tsv").open() as handle:
            trace = list(csv.DictReader(handle, delimiter="\t"))
        energy = float(re.search(r"Final total energy\s*:\s*([-\d.]+)", text)[1])
        hvp = sum(int(row["orbital_hvp_actions"]) for row in trace)
        wall = float(re.search(r"SCF iteration wall time\s*:\s*([\d.]+)", text)[1])
        time = (run / "resource.time").read_text()
        rss = float(re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", time)[1]) / 1024
        previous = old[label]
        rows.append(dict(run=artifact.name, space=label,
                         node=(artifact / "node.txt").read_text().strip(),
                         iterations=int(iterations[1]), orbital_hvp=hvp,
                         energy_Eh=energy, scf_seconds=wall, rss_mib=rss,
                         saved_scf_seconds=float(previous["scf_seconds"]),
                         saved_rss_mib=float(previous["rss_mib"]),
                         saved_iterations=int(previous["iterations"]),
                         saved_orbital_hvp=int(previous["orbital_hvp"]),
                         energy_delta_Eh=energy-float(previous["final_energy_Eh"])))
    save(output, rows)


def generators(artifacts, output):
    rows = []
    for artifact in artifacts:
        for n in (12, 14, 16):
            native_path = artifact / f"{n}-factor-generator.out"
            if not native_path.exists():
                continue
            native = values(native_path)
            reference_path = artifact / f"{n}-factor-generator-reference.out"
            reference = values(reference_path) if reference_path.exists() else {}
            checked = int(reference.get("reference_checked",0))
            def number(key):
                return float(reference[key]) if checked else None
            rows.append(dict(
                run=artifact.name, orbitals=n, node=(artifact / "node.txt").read_text().strip(),
                reference_checked=checked,
                generator_seconds=float(native["warm_generator_seconds"]),
                residual_seconds=float(native["warm_residual_seconds"]),
                rss_mib=float(native["native_max_rss_kib"])/1024,
                paired_generator_seconds=number("warm_generator_seconds"),
                paired_residual_seconds=number("warm_residual_seconds"),
                reference_generator_seconds=number("generator_reference_seconds"),
                reference_sigma_seconds=number("sigma_reference_seconds"),
                reference_cached_residual_seconds=number("residual_generator_reference_seconds"),
                reference_residual_total_seconds=number("residual_total_reference_seconds"),
                generator_error=number("generator_error"),
                residual_error=number("residual_error")))
    save(output, rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="kind", required=True)
    kernel = sub.add_parser("kernels")
    kernel.add_argument("output", type=Path)
    kernel.add_argument("artifacts", type=Path, nargs="+")
    generator = sub.add_parser("generators")
    generator.add_argument("output", type=Path)
    generator.add_argument("artifacts", type=Path, nargs="+")
    chain = sub.add_parser("scf")
    chain.add_argument("artifact", type=Path)
    chain.add_argument("baseline", type=Path)
    chain.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.kind == "kernels":
        kernels(args.artifacts, args.output)
    elif args.kind == "generators":
        generators(args.artifacts, args.output)
    else:
        scf(args.artifact, args.baseline, args.output)


if __name__ == "__main__":
    main()
