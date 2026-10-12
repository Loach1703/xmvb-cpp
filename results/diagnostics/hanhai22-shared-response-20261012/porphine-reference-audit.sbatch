#!/usr/bin/env bash
#SBATCH --job-name=porphine-reference-audit
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=porphine-reference-audit-%j.out
#SBATCH --error=porphine-reference-audit-%j.err
set -eo pipefail
. /etc/profile
module purge
module load anaconda/3/22.05
root="$SLURM_SUBMIT_DIR"
artifact="$root/porphine-reference-audit-$SLURM_JOB_ID"
mkdir -p "$artifact/python"
python3 -m zipfile -e "$root/pyscf-2.14.0-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl" "$artifact/python"
export PYTHONPATH="$artifact/python"
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 -u "$root/porphine-reference-audit.py" "$root" "$@" > "$artifact/audit.out" 2> "$artifact/audit.err"
