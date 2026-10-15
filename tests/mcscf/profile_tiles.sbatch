#!/usr/bin/env bash
#SBATCH --job-name=tile-profile
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=00:05:00
#SBATCH --output=tile-profile-%j.out
#SBATCH --error=tile-profile-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
deps=/home/guqqgroup/xiatao/software/xmvb-deps-20261012
mkdir -p "$source_dir/tile-profile-$SLURM_JOB_ID"
cd "$source_dir/tile-profile-$SLURM_JOB_ID"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH="$deps/lib:${LD_LIBRARY_PATH:-}"
hostname > node.txt
lscpu > cpu.txt
sha256sum "$source_dir/build/src/test_structure_tiles" > binary.sha256
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 env XMVB_PHASE_TIMING=1 /usr/bin/time -v -o native.time \
 "$source_dir/build/src/test_structure_tiles" --selected12 132 > native.out 2> native.err
srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
 bash -c 'command -v gprofng' > profiler.txt 2> profiler.err || exit 0
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 gprofng collect app -p hi -o capture.er \
 "$source_dir/build/src/test_structure_tiles" --selected12 132 > capture.out 2> capture.err
srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
 gprofng display text -header -statistics capture.er > quality.txt
srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
 gprofng display text -limit 35 -functions -calltree capture.er > hotspots.txt
