"""Bounded, same-fixture PySCF AH thread diagnosis (one macro iteration)."""

import argparse
from collections import Counter
import ctypes.util
import json
import os
from pathlib import Path
import resource
import runpy
import sys
import threading
import time

if not os.environ.get("SLURM_JOB_ID"):
    raise SystemExit("Slurm compute allocation required")

affinity_before_numpy = sorted(os.sched_getaffinity(0))
import numpy as np
import scipy
from threadpoolctl import threadpool_info, threadpool_limits
from pyscf import lib, mcscf


def libraries():
    names = ("blas", "mkl", "gomp", "iomp", "omp", "pyscf", "np_helper", "libcint")
    paths = set()
    for line in Path("/proc/self/maps").read_text().splitlines():
        candidate = line.split()[-1]
        if "/" in candidate and any(name in candidate.lower() for name in names):
            paths.add(candidate)
    return sorted(paths)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--blas", type=int, choices=(1, 32), required=True)
parser.add_argument("--omp", type=int, choices=(1, 32), default=1)
parser.add_argument("--original", type=Path, required=True)
parser.add_argument("--fixture", type=Path, required=True)
parser.add_argument("--metadata", type=Path, required=True)
parser.add_argument("--source", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()

report = {"blas_requested": args.blas,
          "omp_requested": args.omp,
          "affinity_before_numpy": affinity_before_numpy,
          "affinity": sorted(os.sched_getaffinity(0)),
          "slurm_job_id": os.environ["SLURM_JOB_ID"],
          "numpy_version": np.__version__, "scipy_version": scipy.__version__,
          "find_library": {name: ctypes.util.find_library(name)
                           for name in ("blas", "openblas", "mkl_rt", "gomp", "omp")},
          "threadpool_before": threadpool_info(),
          "pyscf_omp_before": lib.num_threads()}

original_factory = mcscf.CASSCF


def capped_factory(*factory_args, **factory_kwargs):
    mc = original_factory(*factory_args, **factory_kwargs)
    original_kernel = mc.kernel

    def one_macro(*kernel_args, **kernel_kwargs):
        mc.max_cycle_macro = 1
        lib.num_threads(args.omp)
        report["pyscf_omp_during"] = lib.num_threads()
        report["threadpool_prelimit"] = threadpool_info()
        with threadpool_limits(limits=args.blas, user_api="blas"):
            report["threadpool_during"] = threadpool_info()
            report["libraries_during"] = libraries()
            report["rss_before_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            stop = threading.Event()
            patterns = Counter()
            maximum = {"count": 0, "masks": {}}

            def sample_threads():
                while not stop.is_set():
                    masks = Counter()
                    for entry in Path("/proc/self/task").iterdir():
                        try:
                            mask = tuple(sorted(os.sched_getaffinity(int(entry.name))))
                        except ProcessLookupError:
                            continue
                        masks[",".join(map(str, mask))] += 1
                    total = sum(masks.values())
                    patterns[json.dumps(dict(masks), sort_keys=True)] += 1
                    if total > maximum["count"]:
                        maximum["count"] = total
                        maximum["masks"] = dict(masks)
                    stop.wait(0.1)

            sampler = threading.Thread(target=sample_threads, daemon=True)
            sampler.start()
            wall = time.monotonic()
            cpu = time.process_time()
            try:
                return original_kernel(*kernel_args, **kernel_kwargs)
            finally:
                stop.set()
                sampler.join()
                report["worker_affinity_max_sample"] = maximum
                report["worker_affinity_patterns"] = dict(patterns)
                report["kernel_wall_seconds"] = time.monotonic() - wall
                report["kernel_cpu_seconds"] = time.process_time() - cpu
                report["rss_after_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                report["threadpool_after_kernel"] = threadpool_info()

    mc.kernel = one_macro
    return mc


mcscf.CASSCF = capped_factory
sys.argv = [str(args.original), str(args.fixture), str(args.metadata),
            str(args.source), str(args.output), "--solver", "ah", "--expected-case", "FeCl2"]
try:
    runpy.run_path(str(args.original), run_name="__main__")
except RuntimeError as error:
    if not str(error).startswith("PySCF ah CASSCF did not converge"):
        raise
    report["expected_one_macro_nonconvergence"] = True
finally:
    report["libraries_after"] = libraries()
    report["threadpool_after"] = threadpool_info()
    (args.output.parent / (args.output.name + "-diagnostic.json")).write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
