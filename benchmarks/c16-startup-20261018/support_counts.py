#!/usr/bin/env python3
"""Count symbolic support, not amplitudes, for half-filled singlet COV spaces.

This is a combinatorial execution-plan analysis, not an electronic-structure
calculation or a runtime/memory measurement. It mirrors the prefix condition
of ProductPaths' complete-covalent specialization.
"""

import csv
import math
import sys


def support_count(n: int, depth: int) -> int:
    occupied = n // 2
    counts = {(0, 0): 1}
    for site in range(n):
        following = {}
        for (alpha, beta), multiplicity in counts.items():
            for da, db in ((0, 0), (0, 1), (1, 0), (1, 1)):
                a, b = alpha + da, beta + db
                if a <= occupied and b <= occupied and a + b >= site + 1 - depth:
                    following[a, b] = following.get((a, b), 0) + multiplicity
        counts = following
    return counts[occupied, occupied]


def main() -> None:
    writer = csv.writer(sys.stdout)
    writer.writerow(("active", "structures", "strings_per_spin", "products",
                     "support", "one_body_halo", "dense_vector_MiB",
                     "compact_input_MiB", "accepted_two_images_MiB"))
    for n in map(int, sys.argv[1:] or (12, 14, 16)):
        if n < 2 or n % 2:
            raise ValueError("expected a positive, even, half-filled orbital count")
        strings = math.comb(n, n // 2)
        support = support_count(n, 0)
        writer.writerow((n, strings // (n // 2 + 1), strings, strings**2,
                         support, support_count(n, 1),
                         f"{8 * strings**2 / 2**20:.6f}",
                         f"{12 * support / 2**20:.6f}",
                         f"{20 * support / 2**20:.6f}"))


if __name__ == "__main__":
    main()
