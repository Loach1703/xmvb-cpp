#!/usr/bin/env bash
#SBATCH --job-name=simd-kernel
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15-00:00:00
#SBATCH --output=simd-kernel-%j.out
#SBATCH --error=simd-kernel-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
baseline="$1"
updated="$2"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export PATH="$dependencies/bin:$PATH"
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
root="$SLURM_SUBMIT_DIR/simd-kernel-$SLURM_JOB_ID"
mkdir -p "$root"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
for variant in baseline updated; do
    source_dir="$baseline"
    [[ "$variant" != updated ]] || source_dir="$updated"
    # Identical benchmark harness, original library sources in baseline.
    cmake --build "$source_dir/build" --target test_orthogonal_direct_ci -j16 \
        > "$root/$variant-build.log" 2>&1
    exe="$source_dir/build/src/test_orthogonal_direct_ci"
    sha256sum "$exe" > "$root/$variant-hash.txt"
    for threads in 1 8 32; do
        srun --exact --cpus-per-task="$threads" --cpu-bind=cores \
            env OMP_NUM_THREADS="$threads" /usr/bin/time -v \
            -o "$root/$variant-$threads.time" "$exe" --benchmark \
            > "$root/$variant-$threads.csv" 2> "$root/$variant-$threads.err"
    done
done
