#!/usr/bin/env python3
"""Compare saved fixed-orbital probes; never executes a numerical binary."""

import argparse
import csv
from pathlib import Path


def fields(path):
    return dict(line.split("=", 1) for line in path.read_text().splitlines() if "=" in line)


def save(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("saved", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    history=[]
    for name in ("source_response.csv","source_state_response.csv","link_layout_responses.csv"):
        with (args.saved / name).open() as stream:
            history.extend(csv.DictReader(stream))
    with (args.saved / "source_cover.csv").open() as stream:
        energies = {int(row["case"].split("-")[-1]): float(row["energy_hartree"])
                    for row in csv.DictReader(stream) if row["case"].startswith("selected12-")}
    rows = []
    node = (args.run / "node.txt").read_text().strip()
    for count in (1, 3, 20, 132):
        jet = fields(args.run / f"{count}-product-jet12.out")
        response = fields(args.run / f"{count}-structure-response12.out")
        action = fields(args.run / f"{count}-selected12.out")
        old = [row for row in history if row["orbitals"] == "12" and int(row["structures"]) == count]
        best = min(old, key=lambda row: float(row["warm_seconds"]))
        latest = max(old,key=lambda row:int(row["run"].rsplit("-",1)[-1]))
        new_time = float(response["warm_response_seconds"])
        rows.append(dict(
            run=args.run.name, node=node, orbitals=12, structures=count,
            scope="fixed_orbital_action_and_derivatives_not_converged_SCF",
            rayleigh_energy_Eh=float(action["rayleigh_energy"]),
            saved_energy_delta_Eh=float(action["rayleigh_energy"])-energies[count],
            accepted_s=float(jet["accepted_seconds"]),
            orbital_direction_s=float(jet["direction_seconds"]),
            cached_reverse_s=float(jet["hvp_seconds"]),
            hvp_fd_error=float(jet["hvp_fd_error"]),
            gradient_fd_error=float(jet["lagrangian_fd_error"]),
            accepted_image_cache_bytes=int(jet["accepted_image_cache_bytes"]),
            jet_probe_rss_mib=int(jet["native_max_rss_kib"])/1024,
            response_s=new_time, response_rss_mib=int(response["native_max_rss_kib"])/1024,
            response_checksum=float(response["checksum"]),
            saved_checksum_delta=float(response["checksum"])-float(latest["checksum"]),
            latest_saved_run=latest["run"], latest_saved_response_s=float(latest["warm_seconds"]),
            latest_saved_response_rss_mib=float(latest["rss_mib"]),
            best_saved_run=best["run"], best_saved_response_s=float(best["warm_seconds"]),
            best_saved_response_rss_mib=float(best["rss_mib"]),
            speedup_over_best_saved=float(best["warm_seconds"])/new_time,
            hs_warm_s=float(action["warm_apply_seconds"]),
            hs_rss_mib=int(action["max_rss_kib"])/1024))
    save(args.output, rows)
    for row in rows:
        print(row["structures"], "response", row["response_s"], "RSS", row["response_rss_mib"],
              "best_saved", row["best_saved_response_s"], "energy_delta", row["saved_energy_delta_Eh"],
              "HVP_FD", row["hvp_fd_error"])


if __name__ == "__main__":
    main()
