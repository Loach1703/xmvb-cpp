#!/usr/bin/env bash
#SBATCH --job-name=C12-ci-model
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=64G
#SBATCH --output=C12-ci-model-%j.out
#SBATCH --error=C12-ci-model-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
set -u
source_dir="$SLURM_SUBMIT_DIR"
saved=$(realpath -e "$1")
artifact="$source_dir/c12-ci-model-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
cmake --build "$source_dir/build" --target test_selected_metric_point -j16 > "$artifact/build.log" 2>&1
cp "$source_dir/build/src/test_selected_metric_point" "$artifact/"
cp "$saved/cc-pvdz.gbs" "$saved/input.orb" "$artifact/"
# Inspect the saved endpoint's operator, not a rerun of its SCF optimization.
python3 - "$saved/input.xmi" "$artifact/input.orb" "$artifact/input.xmi" <<'PY'
from pathlib import Path
import re, sys
source, orbitals, target=map(Path,sys.argv[1:])
text,count=re.subn(r'(\$GUS\b).*?(\$END\b)',
    lambda m:m[1]+'\n'+orbitals.read_text()+'\n'+m[2],source.read_text(),flags=re.S|re.I)
assert count==1
target.write_text(text)
PY
cd "$artifact"
hostname > node.txt
lscpu > cpu.txt
sha256sum test_selected_metric_point input.xmi input.orb cc-pvdz.gbs > hashes.sha256
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
grep Cpus_allowed_list /proc/self/status > affinity.txt
/usr/bin/time -v -o resource.time ./test_selected_metric_point input.xmi --ci-model > models.out 2> models.err
RUN
