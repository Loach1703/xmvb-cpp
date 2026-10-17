# Rejected primary-thread callback ownership

Job 1264637 on Hanhai22 anode40 used the same 32 cores, original CAS12/LOFLEA
input and unchanged precision. Every callback ran on team thread0, verified
by the phase min/max owner fields. Original barriers were preserved explicitly.
All 9 tests passed and both native/phase calculations converged in 4 steps /
36 HVPs with matching energies.

| Metric | Saved address plan | Fixed primary |
|---|---:|---:|
| Native wall / s | 30.863387 | 32.613109 |
| Phase-enabled wall / s | 31.477191 | 33.024193 |
| Native RSS / MiB | 2824.676 | 2714.836 |
| Phase RSS / MiB | 2631.992 | 2603.371 |

The experiment lowers measured RSS but increases time about 5%. It does not
establish that the RSS difference is specifically malloc-arena retention;
placement and task scheduling also change. Production code was restored to
the existing single-region policy, and the added owner counters were removed.
The candidate patch, hashes and raw logs are evidence only, not an alternate
maintained implementation. No old executable was rerun.
