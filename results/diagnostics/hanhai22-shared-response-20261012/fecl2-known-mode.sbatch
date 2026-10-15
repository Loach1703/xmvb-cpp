#!/usr/bin/env bash
#SBATCH --job-name=fecl2-known-mode
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --array=0-1
#SBATCH --output=fecl2-mode-%A_%a.out
#SBATCH --error=fecl2-mode-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
artifact="$root/fecl2-mode-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$artifact"
cp "$root/fecl2-first-step-1247319/input.xmi" "$root/fecl2-first-step-1247319/def2-tzvp.gbs" "$artifact/"
cd "$artifact"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
override=()
if (( SLURM_ARRAY_TASK_ID == 1 )); then
  override=(--orbital-value-table-bin "$root/fecl2-first-step-1247319/first-step.bin")
fi
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  "$root/fecl2-update-1247340/check_exact_ctx_hvp" input.xmi --eigensolver davidson \
  --packed-direction-bin "$root/known-negative-direction.bin" \
  --response-tolerance 0.00000001 --step 0.0001 "${override[@]}" \
  > output.out 2> run.err
