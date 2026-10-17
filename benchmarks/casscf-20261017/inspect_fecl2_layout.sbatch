#!/usr/bin/env bash
#SBATCH --job-name=fe-layout
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:03:00
#SBATCH --output=fe-layout-%j.out
#SBATCH --error=fe-layout-%j.err
set -eo pipefail
. /etc/profile
module purge
module load anaconda3/25.06
set -u
layout=/home/guqqgroup/taoxia/xmvb-runs/oeo-pyscf-20261011/runs-282255/FeCl2/cpp/input.molden
basis=/home/guqqgroup/taoxia/work/xmvb-ciah-subspace-20261017/basis/def2-tzvp.gbs
wheel=/home/guqqgroup/taoxia/xmvb-runs/oeo-pyscf-20261011/wheels/pyscf-2.14.0-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
overlay="$SLURM_SUBMIT_DIR/pyscf-overlay-$SLURM_JOB_ID"
mkdir "$overlay"
/opt/anaconda3/25.06/bin/python3 -m zipfile -e "$wheel" "$overlay"
export PYTHONPATH="$overlay" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
srun --exact --ntasks=1 --cpus-per-task=1 \
  /opt/anaconda3/25.06/bin/python3 -u "$SLURM_SUBMIT_DIR/inspect_fecl2_layout.py" \
  "$layout" "$basis"
