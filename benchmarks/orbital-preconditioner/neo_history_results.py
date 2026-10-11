"""Extract historical-strategy evidence from saved native NEO runs."""

import argparse
import csv
import re
from pathlib import Path


def extract(directory):
    output = (directory / "output.out").read_text()
    resource = (directory / "resource.time").read_text()
    trace = list(csv.DictReader((directory / "trace.tsv").open(), delimiter="\t"))

    def field(pattern, text):
        match = re.search(pattern, text)
        if not match:
            raise ValueError(f"missing {pattern} in {directory}")
        return match.group(1)

    if "Exit status: 0" not in resource:
        raise ValueError(f"failed run: {directory}")
    hashes = (directory / "hashes.txt").read_text().splitlines()
    return dict(
        run=str(directory),
        reported_steps=int(field(r"VBSCF converged in\s+(\d+) iterations", output)),
        accepted_points=sum(int(row["keyframes"]) for row in trace),
        hvp=sum(int(row["orbital_hvp_actions"]) for row in trace),
        structure_actions=sum(int(row["structure_response_actions"]) for row in trace),
        energy=float(field(r"Final total energy\s*:\s*([-\d.]+)", output)),
        physical_gradient=float(trace[-1].get("accepted_physical_gradient_norm", "nan")),
        wall_s=float(field(r"End-to-end wall time\s*:\s*([\d.]+)", output)),
        rss_mib=int(field(r"Maximum resident set size \(kbytes\):\s*(\d+)", resource)) / 1024,
        executable_sha256=hashes[0].split()[0],
        input_sha256=hashes[1].split()[0],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("roots", nargs="+", type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    args = parser.parse_args()
    rows = [extract(output.parent) for root in args.roots
            for output in sorted(root.rglob("output.out"))]
    if not rows:
        raise ValueError("no native runs found")
    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
