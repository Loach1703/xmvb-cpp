#!/usr/bin/env python3
"""Parse opt-in matrix shapes; counts are operation sites, not intercepted library calls."""

import argparse
import csv
from pathlib import Path

BACKENDS = {
    "neo_metric_factor": "Eigen LLT / LAPACKE dpotrf",
    "neo_structure_factor": "Eigen SelfAdjointEigenSolver / LAPACKE dsyev (nontrivial sizes)",
    "neo_spectrum": "explicit LAPACKE dsyevd",
    "chart_constraint_svd": "Eigen JacobiSVD / LAPACKE dgesvd under current build macros",
    "chart_gauge_svd": "Eigen JacobiSVD / LAPACKE dgesvd under current build macros",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for directory in sorted(args.root.iterdir()):
        log = directory / "run.err"
        resources = directory / "resource.time"
        if not log.is_file() or not resources.is_file():
            continue
        if "Exit status: 0" not in resources.read_text():
            continue
        for line in log.read_text().splitlines():
            fields = line.split(",")
            if len(fields) != 7 or fields[0] != "XMVB_MATRIX" or fields[1] == "stage":
                continue
            count, min_rows, max_rows, min_cols, max_cols = map(int, fields[2:])
            rows.append(dict(
                variant=directory.name,
                stage=fields[1],
                backend=BACKENDS.get(fields[1], "not mapped"),
                operation_calls=count,
                min_rows=min_rows, max_rows=max_rows,
                min_cols=min_cols, max_cols=max_cols,
            ))
    if not rows:
        raise RuntimeError("no completed matrix-shape records found")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {len(rows)} operation-shape records to {args.output}")
    print("Backend mapping is source-verified for the current build, not a runtime call interceptor.")
    print("Independent row/column extrema do not specify a joint shape or total workspace size.")


if __name__ == "__main__":
    main()
