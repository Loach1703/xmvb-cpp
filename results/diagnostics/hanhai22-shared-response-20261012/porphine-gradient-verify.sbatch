#!/usr/bin/env bash
#SBATCH --job-name=porphine-gradient-verify
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=00:30:00
#SBATCH --output=porphine-gradient-verify-%j.out
#SBATCH --error=porphine-gradient-verify-%j.err
set -eo pipefail
. /etc/profile
module purge
module load anaconda/3/22.05
root="$SLURM_SUBMIT_DIR"
export PYTHONPATH="$root/porphine-reference-audit-1247977/python" PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=32 OMP_PROC_BIND=false OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
unset OMP_PLACES GOMP_CPU_AFFINITY
run="$root/porphine-gradient-verify-$SLURM_JOB_ID"
mkdir -p "$run"
cd "$run"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time python3 -u "$root/porphine-gradient-verify.py" \
  "$root" "$root/shell-screening-fix-1253607/gradient" > output.out 2> run.err
