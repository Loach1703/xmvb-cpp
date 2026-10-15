#!/usr/bin/env bash
#SBATCH --job-name=tiles-kernel
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=tiles-kernel-%j.out
#SBATCH --error=tiles-kernel-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
deps=/home/guqqgroup/xiatao/software/xmvb-deps-20261012
export LD_LIBRARY_PATH="$deps/lib:${LD_LIBRARY_PATH:-}"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 XMVB_PHASE_TIMING=1
artifact="$SLURM_SUBMIT_DIR/tiles-kernel-$SLURM_JOB_ID"
mkdir -p "$artifact"
hostname > "$artifact/node.txt"
for threads in 8 32; do
  export OMP_NUM_THREADS="$threads"
  srun --exact --ntasks=1 --cpus-per-task="$threads" --cpu-bind=cores /usr/bin/time -v -o "$artifact/$threads.time" \
    "$SLURM_SUBMIT_DIR/build/src/test_structure_tiles" --benchmark12-no-reference \
    > "$artifact/$threads.out" 2> "$artifact/$threads.err"
done
