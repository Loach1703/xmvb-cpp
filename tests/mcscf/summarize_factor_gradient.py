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
            check_path = artifact / f"{n}-direction-check.out"
            directional = values(check_path) if check_path.exists() else {}
            directional_checked = int(directional.get("directional_checked",0))
            phase_path = artifact / f"{n}-direction-check.err"
            plans = re.findall(r"XMVB_FACTOR_GENERATOR,hamiltonian=1,([^\n]+)",
                               phase_path.read_text() if phase_path.exists() else "")
            plan = dict(item.split("=",1) for item in plans[-1].split(",")) if plans else {}
            def number(key):
                return float(reference[key]) if checked else None
            rows.append(dict(
                run=artifact.name, orbitals=n, node=(artifact / "node.txt").read_text().strip(),
                reference_checked=checked,
                directional_checked=directional_checked,
                directional_generator_error=float(directional["directional_generator_error"]) if directional_checked else None,
                directional_residual_error=float(directional["directional_residual_error"]) if directional_checked else None,
                left_tile=plan.get("left_tile"), right_tile=plan.get("right_tile"),
                channel_tile=plan.get("channels"), workspace_bytes=plan.get("workspace_bytes"),
                source_columns=plan.get("source_columns"), source_calls=plan.get("source_calls"),
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


def responses(artifacts, output):
    rows = []
    for artifact in artifacts:
        for path in sorted(artifact.glob("*-structure-response*.out")):
            case = re.fullmatch(r"(\d+)-structure-response(\d+)\.out", path.name)
            if case is None:
                continue
            count, n = map(int, case.groups())
            native = values(path)
            reference_path = artifact / f"{count}-structure-response{n}-reference.out"
            reference = values(reference_path) if reference_path.exists() else {}
            checked = int(reference.get("reference_checked",0))
            log_path = artifact / f"{count}-structure-response{n}.err"
            log = log_path.read_text() if log_path.exists() else ""
            links = re.findall(r"^XMVB_LINKS,.*retained_bytes=(\d+)$",log,re.M)
            build = re.search(r"^XMVB_PHASE,ci_connections,\d+,([\d.]+)$",log,re.M)
            tiles = re.findall(r"^XMVB_TILE,([^\n]+)$",log,re.M)
            tile = dict(item.split("=",1) for item in tiles[0].split(",")) if tiles else {}
            rows.append(dict(
                run=artifact.name, orbitals=n, structures=count,
                node=(artifact / "node.txt").read_text().strip(),
                order="source" if int(native["source_columns"]) else "product",
                cold_seconds=float(native["cold_response_seconds"]),
                warm_seconds=float(native["warm_response_seconds"]) if "warm_response_seconds" in native else None,
                single_sweep=int(native.get("single_sweep",0)),
                input_cache_bytes=native.get("input_resident_bytes"),
                output_cache_bytes=native.get("output_resident_bytes"),
                input_panels=native.get("input_panels"),
                input_evictions=native.get("input_evictions"),
                output_evictions=native.get("output_evictions"),
                integral_adjoint_sweeps=native.get("integral_adjoint_sweeps"),
                retained_link_bytes=sum(map(int,links)) if links else None,
                link_build_seconds=float(build[1]) if build else None,
                first_tile_prepare_seconds=tile.get("prepare"),
                first_tile_read_seconds=tile.get("read"),
                first_tile_compute_seconds=tile.get("compute"),
                first_tile_write_seconds=tile.get("write"),
                rss_mib=float(native["native_max_rss_kib"])/1024,
                prepare_seconds=float(native["prepare_seconds"]),
                source_columns=int(native["source_columns"]),
                source_workspace_bound=int(native["source_workspace_bound"]),
                checksum=float(native["checksum"]),
                reference_checked=checked,
                reference_adjoint_only_seconds=float(reference["reference_adjoint_only_seconds"]) if checked else None,
                image_error=float(reference["image_error"]) if checked else None,
                adjoint_error=float(reference["adjoint_error"]) if checked else None))
    rows.sort(key=lambda row: (row["run"], row["orbitals"], row["structures"]))
    save(output, rows)


def states(artifacts, output):
    rows = []
    for artifact in artifacts:
        for path in sorted(artifact.glob("*-source-state.out")):
            native = values(path)
            n = int(native["orbitals"])
            reference_path = artifact / f"{n}-source-state-reference.out"
            reference = values(reference_path) if reference_path.exists() else {}
            checked = int(reference["reference_checked"]) if reference else 0
            direction = int(native["direction_checked"])
            def measured(key):
                return float(reference[key]) if checked else None
            rows.append(dict(
                run=artifact.name, orbitals=n, structures=int(native["structures"]),
                node=(artifact / "node.txt").read_text().strip(),
                accepted_seconds=float(native["accepted_seconds"]),
                direction_checked=direction,
                direction_seconds=float(native["direction_seconds"]) if direction else None,
                hvp_seconds=float(native["hvp_seconds"]) if direction else None,
                combined_seconds=float(native["direction_seconds"])+float(native["hvp_seconds"]) if direction else None,
                rss_mib=float(native["native_max_rss_kib"])/1024,
                lagrangian_fd_error=float(native["lagrangian_fd_error"]),
                accepted_product_doubles=int(native["accepted_product_doubles"]),
                direction_product_doubles=int(native["direction_product_doubles"]),
                cached_generator_doubles=native.get("cached_generator_doubles"),
                direction_source_columns=native.get("direction_source_columns"),
                direction_workspace_bound=native["direction_workspace_bound"] if direction else None,
                reference_checked=checked,
                gradient_error=measured("gradient_error"), image_error=measured("image_error"),
                hvp_error=measured("hvp_error"),
                reference_direction_seconds=measured("reference_direction_seconds"),
                reference_hvp_seconds=measured("reference_hvp_seconds")))
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
    response = sub.add_parser("responses")
    response.add_argument("output", type=Path)
    response.add_argument("artifacts", type=Path, nargs="+")
    state = sub.add_parser("states")
    state.add_argument("output", type=Path)
    state.add_argument("artifacts", type=Path, nargs="+")
    chain = sub.add_parser("scf")
    chain.add_argument("artifact", type=Path)
    chain.add_argument("baseline", type=Path)
    chain.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.kind == "kernels":
        kernels(args.artifacts, args.output)
    elif args.kind == "generators":
        generators(args.artifacts, args.output)
    elif args.kind == "responses":
        responses(args.artifacts, args.output)
    elif args.kind == "states":
        states(args.artifacts, args.output)
    else:
        scf(args.artifact, args.baseline, args.output)


if __name__ == "__main__":
    main()
