"""Extract the observed orbital displacement for a bounded derivative audit."""
import sys
from pathlib import Path
import numpy as np

steps = Path(sys.argv[1]) / "trace/input/steps"
initial = np.fromfile(steps / "step_000000/orbital_value_table_f64.bin", dtype=np.float64)
first = np.fromfile(steps / "step_000001/orbital_value_table_f64.bin", dtype=np.float64)
assert initial.size == first.size == 4488
(first - initial).tofile(sys.argv[2])
