#!/usr/bin/env bash
#SBATCH --job-name=mpi-probe
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --output=mpi-probe-%j.out
#SBATCH --error=mpi-probe-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module avail 2>&1
module load anaconda3/25.06
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
find "$dependencies" /opt -maxdepth 4 \( -name mpicc -o -name mpicxx -o -name libmpi.so \) -print 2>/dev/null || true
srun --mpi=list
