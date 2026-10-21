#!/usr/bin/env bash
#SBATCH --job-name=pyscf-blas-check
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=00:10:00
#SBATCH --output=pyscf-blas-%j.out
#SBATCH --error=pyscf-blas-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# -eq 2 ]] || exit 64
source_dir=$(realpath -e "$1")
fixture_dir=$(realpath -e "$2")
. /etc/profile
module purge
module load gcc/14.3.0 anaconda3/25.06
set -u
envprefix=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
bundle=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle
wheel=/home/guqqgroup/taoxia/xmvb-runs/oeo-pyscf-20261011/wheels/pyscf-2.14.0-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
python_bin=/opt/anaconda3/25.06/bin/python3
script="$source_dir/benchmarks/casscf-threads-20261017/one_macro.py"
original="$source_dir/benchmarks/casscf-20261017/pyscf_newton.py"
artifact="$SLURM_SUBMIT_DIR/pyscf-blas-${SLURM_JOB_ID}"
mkdir "$artifact"
export LD_LIBRARY_PATH="$bundle/lib:$envprefix/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=1 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export PYTHONUNBUFFERED=1
hostname > "$artifact/node.txt"
lscpu > "$artifact/lscpu.txt"
module list > "$artifact/modules.txt" 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > "$artifact/environment.txt"
sha256sum "$script" "$original" "$source_dir/benchmarks/casscf-20261017/fixture_format.py" \
  "$fixture_dir/ao-mo.bin" "$fixture_dir/metadata.json" "$wheel" > "$artifact/inputs.sha256"
"$python_bin" --version > "$artifact/python-version.txt" 2>&1
mkdir "$artifact/python-overlay"
"$python_bin" -m zipfile -e "$wheel" "$artifact/python-overlay"
export PYTHONPATH="$artifact/python-overlay:$source_dir/benchmarks/casscf-20261017"
for blas in 1 32; do
  srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    /usr/bin/time -v -o "$artifact/blas-${blas}.time" \
    "$python_bin" -u "$script" --blas "$blas" --original "$original" \
    --fixture "$fixture_dir/ao-mo.bin" --metadata "$fixture_dir/metadata.json" \
    --source "$source_dir" --output "$artifact/blas-${blas}" \
    > "$artifact/blas-${blas}.out" 2> "$artifact/blas-${blas}.err"
done
