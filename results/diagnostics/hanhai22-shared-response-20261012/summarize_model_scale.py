"""Summarize opt-in finite-model diagnostics without running electronic structure."""
import argparse
import gzip
import json
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    opener = gzip.open if args.log.suffix == ".gz" else open
    groups = {}
    with opener(args.log, "rt") as stream:
        for line in stream:
            fields = line.strip().split(",")
            if len(fields) < 5 or fields[0] != "NEO_DIAG":
                continue
            if fields[3] not in ("ciah_model_curvature_audit", "ciah_model_scale_audit", "ciah_model_block_audit"):
                continue
            values = {k: float(v) for k, v in (f.split("=", 1) for f in fields[4:])}
            key = (int(fields[1]), values["anchor_norm"])
            group = groups.setdefault(key, dict(macro=key[0], anchor_norm=key[1], curvature=[], scaling=[]))
            if fields[3] == "ciah_model_curvature_audit":
                group["curvature"].append(values)
            elif fields[3] == "ciah_model_scale_audit":
                group["scaling"].append(values)
            else:
                group["blocks"] = values
    for group in groups.values():
        rows = group["scaling"]
        for previous, row in zip(rows, rows[1:]):
            for field in ("energy_error", "gradient_error", "local_energy_error", "local_gradient_error"):
                a, b = abs(previous[field]), abs(row[field])
                row[field+"_order"] = math.log(a/b, previous["fraction"]/row["fraction"]) if a and b else None
    result = json.dumps(list(groups.values()), indent=2)
    if args.output:
        args.output.write_text(result+"\n")
    else:
        print(result)


if __name__ == "__main__":
    main()
