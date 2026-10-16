#!/usr/bin/env bash
#SBATCH --job-name=product-parity
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=product-parity-%j.out
#SBATCH --error=product-parity-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -x ${1:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
artifact="$SLURM_SUBMIT_DIR/product-parity-$SLURM_JOB_ID"
mkdir -p "$artifact"
sha256sum "$1" > "$artifact/binary.sha256"
hostname > "$artifact/node.txt"
for count in 1 3 20 132; do
 for mode in selected12-parity structure-response12-reference; do
  # The complete-carrier oracle exists only in this validation executable.
  # These runs qualify accuracy, not native speed or production memory.
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
   "$1" --"$mode" "$count" > "$artifact/$count-$mode.out" 2> "$artifact/$count-$mode.err"
 done
done
