# Coupled-model admission: remaining core-only systems

## Protocol

Job array `281624`, tasks 0 through 4, ran five independent systems in parallel
on Hanhai25. Each task ran its current/coupled pair sequentially on the same
node with compact 32-core affinity, one BLAS thread, exact Libcint integrals
and Davidson. No optimizer iteration limit was added. Gradient and adjacent
energy thresholds remain 1e-3 and 1e-7 Eh. All ten runs explicitly report
convergence and return zero.

The immutable binaries are those from job `281618`:

- current: `71fe190bd423c658e165693b2c925ccc4ccf5cf2498a4054711a947b854d8807`
- coupled: `8f8d2b10b7c500cd1660144bcf4530d5caa47b692a083fb122da1fcd613e17b5`

The only code change is `neo-admission.patch`. No compilation or production
changes were needed. Saved historical input hashes match within every pair;
initial energies agree within 1e-12 Eh. Times are single matched native runs,
not repeated estimates. Differences of a few percent need caution. RSS is
GNU-time native process maximum, not allocated Slurm memory.

## Results

Each cell is **current core-only / coupled**. NEGPIA is the preceding matched
experiment, job `281618`, not a newly repeated sample.

| System | Points | Orbital HVP | Structure actions | Total s | RSS MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| 240 | 8/6 | 12/12 | 0/179 | 2.031/2.080 | 1770.7/1691.2 |
| 67970 | 7/6 | 13/13 | 0/251 | 3.138/3.445 | 2010.9/1980.1 |
| CONYAH | 7/5 | 13/11 | 0/87 | 2.002/1.952 | 1634.7/1637.2 |
| QOCGUM | 7/6 | 15/16 | 0/317 | 9.119/9.637 | 2088.6/2073.1 |
| YEPKAF | 7/7 | 25/31 | 0/332 | 9.257/11.329 | 2283.3/2204.8 |
| NEGPIA | 17/6 | 168/45 | 0/293 | 19.452/8.426 | 3765.5/3726.3 |

Four new cases reduce accepted points; YEPKAF does not. Only CONYAH has a
slightly lower measured total time among the five new cases; the difference
is small. YEPKAF costs 22.4% more time. Coupled memory is not systematically
higher in these measured samples, but this does not establish a large-space
memory bound.

| System | Current final energy, Eh | Coupled final energy, Eh |
| --- | ---: | ---: |
| 240 | -343.446302347232 | -343.446302343374 |
| 67970 | -442.247231402783 | -442.247231404580 |
| CONYAH | -343.417822678080 | -343.417822700631 |
| QOCGUM | -480.094766640467 | -480.094766678862 |
| YEPKAF | -1137.242583650808 | -1137.242584086473 |
| NEGPIA | -436.262505347363 | -436.262507052760 |

The first four new endpoint differences are below 4e-8 Eh. YEPKAF differs
by 4.36e-7 Eh and ends with a substantially smaller projected gradient in
the coupled run. NEGPIA differs by 1.705e-6 Eh. Both variants satisfy the
unchanged stopping criteria, but these endpoint differences mean the timings
are not comparisons at exactly identical final energy accuracy.

## Interpretation

Coupling removal is a genuine model-quality regression for NEGPIA, but
complete coupling is not a universal performance solution. Most other
core-only cases already need only 12 to 25 orbital HVPs; adding hundreds of
structure actions can offset their reduced nonlinear iteration count.

The cost-only admission rule cannot capture this tradeoff. Equally, a blanket
return to historical admission would introduce new regressions. A replacement
must distinguish omitted-coupling model error from expensive over-certification
and should reuse available response information rather than rebuild a complete
response unconditionally. No system-specific exception, loosened tolerance or
new production admission rule is introduced in this experiment.

The next analysis should compare predicted versus actual gradient changes and
energy decreases along these accepted trajectories, especially NEGPIA before
its radius collapse, to define a model-quality diagnostic before implementation.

Data: `results/history/neo_admission_panel_20261011.csv` and the preceding
`neo_admission_20261011.csv`. Raw logs, input/binary hashes, CPU information,
affinity and native resources are retained in
`results/studies/locality-20261011/neo-admission-panel-281624/` and the matching
remote directory under `/home/guqqgroup/taoxia/xmvb-runs/native-20261011/`.
