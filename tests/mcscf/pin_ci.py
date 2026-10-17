#!/usr/bin/env python3
"""Run an eight-core CI probe within one allocated last-level-cache domain."""

import json
import os
from pathlib import Path
import sys


def cpu_list(text):
    result = set()
    for item in text.strip().split(","):
        limits = item.split("-")
        result.update(range(int(limits[0]), int(limits[-1]) + 1))
    return result


def main():
    threads, record, *command = sys.argv[1:]
    threads = int(threads)
    available = set(os.sched_getaffinity(0))
    domains = set()
    for cpu in sorted(available):
        cache = Path(f"/sys/devices/system/cpu/cpu{cpu}/cache/index3")
        domains.add(tuple(sorted(available & cpu_list((cache / "shared_cpu_list").read_text()))))
    first = command[0] == "--first"
    if first:
        command = command[1:]
    eligible = sorted(domain for domain in domains if len(domain) >= threads)
    if first:
        selected = sorted(available)[:threads]
    else:
        if not eligible:
            raise SystemExit("allocation has no complete requested LLC group")
        selected = eligible[0][:threads]
    os.sched_setaffinity(0, selected)
    with open(record, "a") as stream:
        stream.write(json.dumps(dict(available=sorted(available), cpus=selected,
                                     cache_groups=[domain for domain in sorted(domains)
                                                   if set(domain) & set(selected)],
                                     first=first)) + "\n")
    os.execvp(command[0], command)


if __name__ == "__main__":
    main()
