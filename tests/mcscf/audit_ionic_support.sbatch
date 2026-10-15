#!/usr/bin/env bash
#SBATCH --job-name=ionic-proof
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:02:00
#SBATCH --output=ionic-proof-%j.out
#SBATCH --error=ionic-proof-%j.err
set -eo pipefail
. /etc/profile
module load gcc/14.1.0
artifact="$SLURM_SUBMIT_DIR/ionic-proof-$SLURM_JOB_ID"
mkdir "$artifact"
g++ -O3 -std=c++17 -I src tests/mcscf/audit_ionic_support.cpp -o "$artifact/proof"
srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores "$artifact/proof" > "$artifact/result.txt"
