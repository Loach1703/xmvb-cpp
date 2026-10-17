#!/usr/bin/env bash
#SBATCH --job-name=sparse-strings
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=00:20:00
#SBATCH --output=sparse-%A_%a.out
#SBATCH --error=sparse-%A_%a.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0 anaconda/3/22.05
root="$SLURM_SUBMIT_DIR"
sizes=(12 14 16)
n="${sizes[$SLURM_ARRAY_TASK_ID]}"
run="$root/run-$SLURM_ARRAY_JOB_ID/$n"
mkdir -p "$run"
cd "$run"
cp "$root/probe.cpp" "$root/prepare.py" .
export PYTHONPATH=/home/guqqgroup/xiatao/xmvb-runs/shared-response-20261012/porphine-reference-audit-1247977/python
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
hostname > node.txt
lscpu > cpu.txt
printf '%s\n' 'g++ -std=c++20 -O3 -g -fno-omit-frame-pointer probe.cpp -o probe' > build.txt
srun --exact --cpu-bind=cores g++ -std=c++20 -O3 -g -fno-omit-frame-pointer probe.cpp -o probe
srun --exact --cpu-bind=cores python3 prepare.py "$n" > preparation.out 2> preparation.err
sha256sum probe.cpp prepare.py probe orbitals.txt > hashes.txt
srun --exact --cpu-bind=cores /usr/bin/time -v -o resource.time ./probe > output.out 2> run.err
