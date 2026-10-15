#!/usr/bin/env bash
#SBATCH --job-name=porphine-cross-restart
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=00:30:00
#SBATCH --output=porphine-cross-restart-%A_%a.out
#SBATCH --error=porphine-cross-restart-%A_%a.err
set -eo pipefail
. /etc/profile
module purge
module load anaconda/3/22.05
root="$SLURM_SUBMIT_DIR"
modes=(default-ah newton)
mode="${modes[$SLURM_ARRAY_TASK_ID]}"
artifact="$root/porphine-cross-restart-$SLURM_ARRAY_JOB_ID/$mode"
mkdir -p "$artifact"
cd "$artifact"
export PYTHONPATH="$root/porphine-reference-audit-1247977/python" PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=32 OMP_PROC_BIND=false OMP_MAX_ACTIVE_LEVELS=1
unset OMP_PLACES GOMP_CPU_AFFINITY
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
 python3 -u "$root/porphine-cross-restart.py" "$root" "$mode" > output.out 2> run.err
