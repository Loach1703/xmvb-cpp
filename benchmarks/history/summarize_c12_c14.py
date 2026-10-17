#!/usr/bin/env python3
"""Reproduce the C12/C14 convergence CSV from committed native-run archives."""

import csv
import hashlib
import io
from pathlib import Path
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests/mcscf"))
from summarize_scf_step import elapsed_seconds, field

CASES = (
    ("C12", "historical", "cd958ad7", 1264541,
     "c12-tiles-20261017/evidence.tar.gz", "run-1264541", "before-loflea", "output.out", "input.xmi"),
    ("C12", "current", "d54c9095", 1264965,
     "c12-gemm-orbit-20261017/evidence.tar.gz", "run-1264965", "native", "output.out", "input.xmi"),
    ("C14", "historical", "cd958ad7", 1264476,
     "c14-current-20261017/evidence.tar.gz", "fix-1264476", "", "C14.out", "C14.xmi"),
    ("C14", "current", "c45f9942", 1265059,
     "c14-current-20261017/evidence.tar.gz", "run-1265059", "", "C14.out", "C14.xmi"),
)


def summarize(case):
    system, variant, revision, job, archive, root, subdir, output, input_file = case
    with tarfile.open(ROOT / "results/diagnostics" / archive) as data:
        def read(name):
            with data.extractfile(str(Path(root) / name)) as stream:
                return stream.read()

        base = Path(subdir)
        text = read(base / output).decode()
        resource = read(base / "resource.time").decode()
        trace = list(csv.DictReader(io.StringIO(read(base / "ciah.tsv").decode()), delimiter="\t"))
        rss = field(resource, r"Maximum resident set size \(kbytes\):\s*(\d+)", int)
        row = dict(
            system=system, variant=variant, revision=revision, job=job,
            node=read("node.txt").decode().strip(),
            threads=field(text, r"^CPU threads\s*:\s*(\d+)", int),
            steps=len(trace), hvp=sum(int(r["orbital_hvp_actions"]) for r in trace),
            responses=sum(int(r["structure_response_actions"]) for r in trace),
            per_step_hvp=":".join(r["orbital_hvp_actions"] for r in trace),
            initial_Eh=field(text, r"^Initial total energy\s*:\s*([-\d.]+)"),
            final_Eh=field(text, r"^Final total energy\s*:\s*([-\d.]+)"),
            application_s=field(text, r"^End-to-end wall time\s*:\s*([\d.]+)", float),
            scf_s=field(text, r"^SCF iteration wall time\s*:\s*([\d.]+)", float),
            process_s=elapsed_seconds(resource), rss_KiB=rss, rss_GiB=rss / 1048576,
            input_sha256=hashlib.sha256(read(base / input_file)).hexdigest(),
            termination=field(text, r"^Termination reason\s*:\s*(\S+)"),
        )
    if row["termination"] != "ciah_projected_gradient_tolerance" or row["threads"] != 32:
        raise ValueError(f"unqualified convergence record: {system}/{variant}")
    return row


if __name__ == "__main__":
    rows = [summarize(case) for case in CASES]
    for old, new in zip(rows[::2], rows[1::2]):
        if old["input_sha256"] != new["input_sha256"] or old["initial_Eh"] != new["initial_Eh"]:
            raise ValueError(f"input/initial-energy mismatch: {old['system']}")
    writer = csv.DictWriter(sys.stdout, fieldnames=rows[0].keys(), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
