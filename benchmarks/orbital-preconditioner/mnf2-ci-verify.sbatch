#!/usr/bin/env bash
#SBATCH --job-name=mnf2-ci-verify
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=15-00:00:00
#SBATCH --output=mnf2-ci-verify-%j.out
#SBATCH --error=mnf2-ci-verify-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export PATH="$dependencies/bin:$PATH" CONDA_PREFIX="$dependencies"
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export LC_ALL=C
root="$SLURM_SUBMIT_DIR/mnf2-ci-verify-$SLURM_JOB_ID"
mkdir -p "$root"
hostname > "$root/node.txt"
cmake --build "$source_dir/build" --target xmvb test_orthogonal_direct_ci -j16 \
    > "$root/build.log" 2>&1
for threads in 1 8 32; do
    srun --exact --cpus-per-task="$threads" --cpu-bind=cores \
        env OMP_NUM_THREADS="$threads" "$source_dir/build/src/test_orthogonal_direct_ci" \
        > "$root/test-$threads.out" 2>&1
done
sha256sum "$source_dir/build/src/xmvb-cpp.exe" > "$root/hash.txt"
