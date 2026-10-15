#!/usr/bin/env bash
#SBATCH --job-name=fe-guess-audit
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=guess-audit-%j.out
#SBATCH --error=guess-audit-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run_dir="$SLURM_SUBMIT_DIR/guess-audit-$SLURM_JOB_ID"
mkdir -p "$run_dir"
hostname > "$run_dir/node.txt"
srun --cpu-bind=cores /home/guqqgroup/taoxia/conda-envs/pyscf/bin/python \
  "$SLURM_SUBMIT_DIR/audit_guess.py" "$1" > "$run_dir/audit.json"
