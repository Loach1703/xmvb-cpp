#!/usr/bin/env bash
#SBATCH --job-name=integral-tools
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=00:10:00
#SBATCH --output=tools-%j.out
#SBATCH --error=tools-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0 cmake/3.31.9
set -u
source_dir=$(realpath -e "$1")
artifact="$SLURM_SUBMIT_DIR/tools-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LD_LIBRARY_PATH=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 OMP_MAX_ACTIVE_LEVELS=1
cmake --build "$source_dir/build-final" --target \
  check_active_space_gradient check_exact_two_electron_hvp \
  check_structure_overlap_against_raw_reference check_structure_builder_fast_path \
  check_structure_expansion_signs benchmark_union_graph_single_step \
  validate_multi_leaf_open_state_one_electron -j"$SLURM_CPUS_PER_TASK" > "$artifact/build.log" 2>&1
cd "$artifact"
cp "$source_dir/basis/"{cc-pvdz.gbs,def2-tzvp.gbs} .
for component in overlap one_electron two_electron; do
  "$source_dir/build-final/src/check_active_space_gradient" \
    "$source_dir/testdata/vbscf/F2.xmi" --component "$component" --count 3 \
    > "$component.out" 2> "$component.err"
done
"$source_dir/build-final/src/check_exact_two_electron_hvp" \
  "$source_dir/testdata/vbscf/F2.xmi" --step 1e-5 > exact-jet.out 2> exact-jet.err
python3 - <<'PY'
import re
from pathlib import Path
for component in ("overlap", "one_electron", "two_electron"):
    values = [float(v) for v in re.findall(
        r"abs_error=([\deE+.-]+)", Path(component + ".out").read_text())]
    if len(values) != 3 or max(values) > 1e-6:
        raise RuntimeError(f"{component} finite difference failed: {values}")
    print(component, "maximum absolute error", max(values))
PY
