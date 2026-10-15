#!/usr/bin/env bash
#SBATCH --job-name=porphine-path-analysis
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=porphine-path-analysis-%j.out
#SBATCH --error=porphine-path-analysis-%j.err
set -eo pipefail
. /etc/profile
module purge
module load anaconda/3/22.05
root="$SLURM_SUBMIT_DIR"
run="$root/porphine-path-${1:?Specify completed path array ID}"
args=("$root" "$run")
if [[ $# == 2 ]]; then args+=("$root/porphine-path-$2/pyscf-newton-1"); fi
export PYTHONPATH="$root/porphine-reference-audit-1247977/python" PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
srun --exact --ntasks=1 --cpus-per-task=8 --cpu-bind=cores \
  python3 "$root/porphine-path-analysis.py" "${args[@]}" > "$run/analysis.out" 2> "$run/analysis.err"
