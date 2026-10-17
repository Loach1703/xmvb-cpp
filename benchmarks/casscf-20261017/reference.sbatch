#!/usr/bin/env bash
#SBATCH --job-name=casscf-ah-ref
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --output=casscf-ah-%j.out
#SBATCH --error=casscf-ah-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# -eq 3 ]] || exit 64
source_dir=$(realpath -e "$1")
fixture_dir=$(realpath -e "$2")
case_name=$3
[[ $case_name == F2 || $case_name == FeCl2 ]] || exit 64
[[ -r "$source_dir/parent-revision.txt" && -r "$source_dir/source.sha256" &&
   -r "$fixture_dir/ao-mo.bin" && -r "$fixture_dir/metadata.json" ]] || exit 66
. /etc/profile
module purge
module load gcc/14.3.0 anaconda3/25.06
set -u
envprefix=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
bundle=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle
wheel=/home/guqqgroup/taoxia/xmvb-runs/oeo-pyscf-20261011/wheels/pyscf-2.14.0-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
python_bin=/opt/anaconda3/25.06/bin/python3
script="$source_dir/benchmarks/casscf-20261017/pyscf_newton.py"
artifact="$SLURM_SUBMIT_DIR/casscf-ah-${case_name}-${SLURM_JOB_ID}"
mkdir "$artifact"
export LD_LIBRARY_PATH="$bundle/lib:$envprefix/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export PYTHONUNBUFFERED=1
hostname > "$artifact/node.txt"
lscpu > "$artifact/lscpu.txt"
module list > "$artifact/modules.txt" 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > "$artifact/environment.txt"
cp "$source_dir/parent-revision.txt" "$source_dir/source.sha256" "$artifact/"
(cd "$source_dir" && sha256sum --check source.sha256) \
  > "$artifact/source-verify.out" 2> "$artifact/source-verify.err"
"$python_bin" --version > "$artifact/python-version.txt" 2>&1
mkdir "$artifact/python-overlay"
"$python_bin" -m zipfile -e "$wheel" "$artifact/python-overlay"
export PYTHONPATH="$artifact/python-overlay"
"$python_bin" -c 'import numpy, scipy, h5py, pyscf; print(numpy.__version__, scipy.__version__, h5py.__version__, pyscf.__version__)' \
  > "$artifact/pyscf-version.txt"
sha256sum "$script" "$source_dir/benchmarks/casscf-20261017/fixture_format.py" \
  "$fixture_dir/ao-mo.bin" "$fixture_dir/metadata.json" "$wheel" \
  > "$artifact/input.sha256"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  bash -c 'grep Cpus_allowed_list /proc/self/status' > "$artifact/affinity.txt"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o "$artifact/pyscf-ah.time" \
  "$python_bin" -u "$script" \
  "$fixture_dir/ao-mo.bin" "$fixture_dir/metadata.json" \
  "$source_dir" "$artifact/pyscf-ah" --solver ah --expected-case "$case_name" \
  > "$artifact/pyscf-ah.out" 2> "$artifact/pyscf-ah.err"
