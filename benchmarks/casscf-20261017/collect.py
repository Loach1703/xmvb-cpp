"""Collect completed CASSCF logs; molecular workloads are never rerun."""

import argparse
import csv
import json
import re
from pathlib import Path


def resources(path):
    text = path.read_text() if path.exists() else ""
    wall = re.search(r"Elapsed \(wall clock\).*?: (\S+)", text)
    seconds = 0.0
    if wall:
        for part in wall[1].split(":"):
            seconds = 60 * seconds + float(part)
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    return dict(process_seconds=seconds if wall else "",
                peak_rss_mib=int(rss[1]) / 1024 if rss else "")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = []
    for artifact in sorted(args.root.iterdir()):
        if not artifact.is_dir() or not artifact.name.startswith(("casscf-", "native-", "fixture-")):
            continue
        common = dict(artifact=artifact.name)
        node = artifact / "node.txt"
        common["node"] = node.read_text().strip() if node.exists() else ""
        for result in sorted(artifact.glob("pyscf*/result.json")):
            data = json.loads(result.read_text())
            if not data["converged"]:
                continue
            macros = int(data["trajectory"][-1]["imacro"])
            rows.append(dict(common, case=data["case"], interface="shared_fixture",
                solver=data["solver"], initial_energy=data["initial_casci_energy"],
                energy=data["energy"], macros=macros,
                hvp=data["hvp_actions"], solver_seconds=data["elapsed_seconds"],
                **resources(artifact / (result.parent.name + ".time"))))
        cpp = artifact / "cpp.out"
        if artifact.name.startswith("fixture-"):
            cpp = artifact / "output.out"
        if cpp.exists():
            values = dict(line.split("=", 1) for line in cpp.read_text().splitlines()
                          if "=" in line)
            if values.get("converged") == "1":
                metadata = artifact / "fixture/metadata.json"
                case = (json.loads(metadata.read_text())["case"] if metadata.exists()
                        else artifact.name.split("-")[1])
                time_file = artifact / ("resource.time" if artifact.name.startswith("fixture-")
                                        else "cpp.time")
                rows.append(dict(common, case=case, interface="shared_fixture",
                    solver="C++ coupled CIAH", initial_energy=values["initial_energy"],
                    energy=values["energy"], macros=values["macros"],
                    hvp=values["hvp_actions"], solver_seconds=values["solver_wall_seconds"],
                    **resources(time_file)))
        output = artifact / "output.out"
        if output.exists() and artifact.name.startswith("native-"):
            text = output.read_text()
            if "Converged: yes" not in text:
                continue
            def field(name):
                return re.search(r"^" + re.escape(name) + r": ([^\s]+)",
                                 text, re.MULTILINE)[1]
            rows.append(dict(common, case=artifact.name.split("-")[1], interface="xmi",
                solver="C++ coupled CIAH", initial_energy=field("Initial energy"),
                energy=field("Final energy"), macros=field("Accepted iterations"),
                hvp=field("Hessian-vector actions"), solver_seconds=field("Solve time"),
                **resources(artifact / "resource.time")))
    fields = ("case", "interface", "solver", "initial_energy", "energy", "macros",
              "hvp", "solver_seconds", "process_seconds", "peak_rss_mib", "node", "artifact")
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Collected {len(rows)} completed results into {args.output}")


if __name__ == "__main__":
    main()
