# Fe porphine CIAH status, 2026-10-12 approximately 19:26 (UTC+08)

| Quantity | Singlet | Triplet |
|---|---:|---:|
| State | Running | Failed, exit 2 |
| Accepted outer steps | 46 | 28 |
| Recorded accepted-step orbital HVP subtotal | 466 | 331 |
| Latest/final energy (Eh) | -2244.741039162000 | -2244.790784197357 |
| Existing PySCF final energy (Eh) | -2244.734377710447 | -2244.782671264980 |
| CPP minus PySCF (mEh) | -6.661452 | -8.112932 |
| Elapsed wall time | About 56 min, ongoing | 2375.82 s |
| Peak RSS (GiB) | 38.620 observed so far | 38.632 final |

Triplet terminated with `ciah_macro_energy_rejected`, not a timeout or OOM.
It did not satisfy the normal optimizer stopping criteria. Singlet remains
active. Neither result establishes agreement with the existing PySCF endpoint.
Being lower in energy does not by itself establish which endpoint is correct.

HVP counts are sums of the accepted outer-iteration trace rows. They exclude
actions during an unfinished or failed macrostep and are not certified total
action counts for these runs.

Triplet wall/RSS are from native GNU `time -v`: peak 40,508,472 KiB.
Singlet peak is the running process's `/proc/PID/status` `VmHWM`, sampled at
19:25:23 by read-only Slurm job **1247902** on anode56: 40,496,196 KiB.
It is a peak observed so far, not the final peak. The probe confirmed 32 threads
and matched the process working directory to this exact singlet job.

The probe was submitted separately because login-node `sstat` and an overlapping
`srun` both failed to resolve anode56 through the local Slurm configuration.
The original chemistry jobs were not restarted or changed. Raw output and
traces are mirrored under `porphine-ciah-1247745/`.
