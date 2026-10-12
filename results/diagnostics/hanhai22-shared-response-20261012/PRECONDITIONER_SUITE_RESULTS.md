# Suite 1247682: partial results

Snapshot on 2026-10-12: 29 of 30 runs finished; LOFLEA one-electron remains
running. Thirteen paired systems have both variants nominally converged.
CERRAS has only analytic-diagonal convergence; FeCl2 one-electron stops on the
known higher-energy branch. No numerical tolerances were changed.

| System | Preconditioner | Status | Steps | HVP | Energy (Eh) | Wall (s) | RSS (MiB) |
|---|---|---|---:|---:|---:|---:|---:|
| F2 | hessian-diagonal | converged | 3 | 24 | -198.751155830525 | 0.28 | 17.1 |
| F2 | one-electron | converged | 3 | 23 | -198.751155830526 | 0.27 | 17.2 |
| 241 | hessian-diagonal | converged | 3 | 26 | -230.720590392852 | 2.63 | 798.0 |
| 241 | one-electron | converged | 3 | 23 | -230.720590392889 | 2.58 | 875.3 |
| 240 | hessian-diagonal | converged | 3 | 40 | -343.446302353192 | 5.90 | 1550.0 |
| 240 | one-electron | converged | 4 | 37 | -343.446302353107 | 6.36 | 1536.6 |
| MnF2 | hessian-diagonal | converged | 11 | 169 | -1348.893353223712 | 5.98 | 463.4 |
| MnF2 | one-electron | converged | 18 | 288 | -1348.893353118897 | 9.91 | 472.2 |
| FeCl2 | hessian-diagonal | converged | 21 | 453 | -2181.624464717314 | 13.30 | 398.3 |
| FeCl2 | one-electron | converged | 22 | 550 | -2181.617644130828 | 14.08 | 418.0 |
| 7975 | hessian-diagonal | converged | 3 | 29 | -285.733616934357 | 4.05 | 1307.4 |
| 7975 | one-electron | converged | 3 | 27 | -285.733616934361 | 4.03 | 1433.2 |
| 7963 | hessian-diagonal | converged | 3 | 28 | -285.731743368988 | 4.47 | 1434.2 |
| 7963 | one-electron | converged | 3 | 32 | -285.731743369378 | 4.53 | 1456.0 |
| YAMSAI | hessian-diagonal | converged | 3 | 37 | -305.561919470877 | 2.60 | 780.9 |
| YAMSAI | one-electron | converged | 3 | 30 | -305.561919471035 | 2.48 | 797.3 |
| NEGPIA | hessian-diagonal | converged | 4 | 73 | -436.262507343068 | 19.07 | 3570.4 |
| NEGPIA | one-electron | converged | 6 | 112 | -436.262507342621 | 25.17 | 3588.4 |
| QOCGUM | hessian-diagonal | converged | 4 | 55 | -480.094766681825 | 15.00 | 1926.1 |
| QOCGUM | one-electron | converged | 4 | 56 | -480.094766681838 | 16.65 | 1827.5 |
| 67970 | hessian-diagonal | converged | 4 | 53 | -442.247231416377 | 10.12 | 1967.2 |
| 67970 | one-electron | converged | 4 | 46 | -442.247231415439 | 10.07 | 1951.3 |
| CERRAS | hessian-diagonal | converged | 5 | 71 | -397.093435003879 | 42.33 | 1836.0 |
| CERRAS | one-electron | not converged | 4 | 76+ | -397.093435016820 | 61.59 | 1850.6 |
| LOFLEA | hessian-diagonal | converged | 4 | 46 | -399.106646007257 | 94.99 | 2735.3 |
| LOFLEA | one-electron | running | — | — | — | — | — |
| 10698 | hessian-diagonal | converged | 4 | 53 | -422.611824623414 | 20.60 | 5092.6 |
| 10698 | one-electron | converged | 4 | 47 | -422.611824623179 | 20.39 | 4972.7 |
| F2-OEO | hessian-diagonal | converged | 4 | 18 | -198.689800055622 | 0.25 | 17.6 |
| F2-OEO | one-electron | converged | 4 | 29 | -198.689800054547 | 0.36 | 17.9 |

On failed runs, steps are accepted steps; HVP marked '+' is only the recorded subtotal.

Paired energy differences (one-electron minus analytic diagonal):
- F2: -9.94759830064e-13 Eh
- 241: -3.70050656784e-11 Eh
- 240: +8.49809111969e-11 Eh
- MnF2: +1.04814944279e-07 Eh
- FeCl2: +0.0068205864859 Eh
- 7975: -4.03588273912e-12 Eh
- 7963: -3.90002696804e-10 Eh
- YAMSAI: -1.58024704433e-10 Eh
- NEGPIA: +4.47016645921e-10 Eh
- QOCGUM: -1.30171429191e-11 Eh
- 67970: +9.38030098041e-10 Eh
- 10698: +2.34990693571e-10 Eh
- F2-OEO: +1.07499431579e-09 Eh
