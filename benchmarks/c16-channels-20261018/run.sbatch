#!/usr/bin/env bash
#SBATCH --job-name=C16-channels
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=channels-%j.out
#SBATCH --error=channels-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0
set -u
source_dir=$(realpath -e "$1")
run_dir="$SLURM_SUBMIT_DIR/run-$SLURM_JOB_ID"
mkdir -p "$run_dir"
export OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
g++ -std=c++23 -O3 -DNDEBUG -g1 -fopenmp -I"$source_dir/src" \
  "$source_dir/benchmarks/c16-channels-20261018/audit.cpp" -o "$run_dir/audit"
hostname > "$run_dir/node.txt"
lscpu > "$run_dir/cpu.txt"
sha256sum "$run_dir/audit" "$source_dir/benchmarks/c16-channels-20261018/audit.cpp" > "$run_dir/hashes.txt"
# Exhaustive small-space envelope checks, then a deterministic C16 census.
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  "$run_dir/audit" 8 1 > "$run_dir/c8.csv"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o "$run_dir/c16.time" "$run_dir/audit" 16 257 > "$run_dir/c16.csv"
