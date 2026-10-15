#!/usr/bin/env bash
#SBATCH --job-name=pyrrole-pyscf
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=pyscf-%j.out
#SBATCH --error=pyscf-%j.err
set -eo pipefail
. /etc/profile
module purge
module load anaconda/3/22.05
root="$SLURM_SUBMIT_DIR"
mkdir -p "$root/pyscf-$SLURM_JOB_ID"
cd "$root/pyscf-$SLURM_JOB_ID"
cp "$root/inputs/oeo.xmi" input.xmi
cp "$root/run-1255187/oeo-block-lbfgs/input.molden" basis_layout.molden
cp "$root/pyscf_optimize.py" .
export PYTHONPATH=/home/guqqgroup/xiatao/xmvb-runs/shared-response-20261012/porphine-reference-audit-1247977/python
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=32 OMP_PROC_BIND=false
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
unset OMP_PLACES GOMP_CPU_AFFINITY
hostname > node.txt
sha256sum input.xmi basis_layout.molden pyscf_optimize.py > hashes.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time python3 -u pyscf_optimize.py > output.out 2>run.err || status=$?
printf '%s\n' "$status" >run.exit
exit "$status"
