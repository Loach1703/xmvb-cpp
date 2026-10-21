#!/usr/bin/env python3
"""Join current exact/RI measurements, requiring the same binary and CPU mask."""

import argparse
import csv
from pathlib import Path


def read(path):
    with path.open() as stream:
        return {row["case"]: row for row in csv.DictReader(stream)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ri", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("exact", type=Path, nargs="+")
    args = parser.parse_args()
    exact = {}
    for path in args.exact:
        exact.update(read(path))
    rows = []
    keys = ("job", "steps", "hvp", "ci_actions", "process_s", "rss_GiB",
            "hvp_s", "effective_cores", "initial_Eh", "final_Eh")
    for case, ri in read(args.ri).items():
        row = dict(case=case, status="missing_current_exact",
                   node=ri["node"], affinity=ri["affinity"],
                   binary_sha256=ri["binary_sha256"], input_sha256=ri["input_sha256"])
        for name in ("ri", "exact"):
            run = ri if name == "ri" else exact.get(case, {})
            row.update({f"{name}_{key}": run.get(key, "") for key in keys})
            row[f"{name}_hvp_mean_s"] = float(run["hvp_s"])/int(run["hvp"]) if run else ""
        row["exact_over_ri_time"] = row["ri_minus_exact_mEh"] = ""
        if case in exact:
            reference = exact[case]
            for key in ("node", "affinity", "binary_sha256", "input_sha256", "basis"):
                if ri[key] != reference[key]:
                    raise RuntimeError(f"{case}: unmatched {key}")
            if ri["converged"] != "True" or reference["converged"] != "True":
                raise RuntimeError(f"{case}: unconverged exact/RI comparison")
            row["status"] = "matched"
            row["exact_over_ri_time"] = float(reference["process_s"])/float(ri["process_s"])
            row["ri_minus_exact_mEh"] = 1000*(float(ri["final_Eh"])-float(reference["final_Eh"]))
        rows.append(row)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
