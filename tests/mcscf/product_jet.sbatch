#!/usr/bin/env bash
#SBATCH --job-name=product-jet
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:20:00
#SBATCH --output=product-jet-%j.out
#SBATCH --error=product-jet-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
artifact="$source_dir/product-jet-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
env | sort > "$artifact/environment.txt"
cmake -S "$source_dir" -B "$source_dir/build" > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build" --target test_structure_tiles -j8 > "$artifact/build.log" 2>&1
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
cp "$source_dir/build/src/test_structure_tiles" "$artifact/"
sha256sum "$artifact/test_structure_tiles" > "$artifact/binary.sha256"
find "$source_dir/src/mcscf/ci" "$source_dir/src/vbscf/structures/assembly" -type f -print0 |
 sort -z | xargs -0 sha256sum > "$artifact/sources.sha256"
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 bash -c 'grep Cpus_allowed_list /proc/self/status' > "$artifact/affinity.txt"
for count in 1 3 20 132; do
 for mode in product-jet12 structure-response12 selected12; do
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
   /usr/bin/time -v -o "$artifact/$count-$mode.time" \
   "$artifact/test_structure_tiles" --"$mode" "$count" \
   > "$artifact/$count-$mode.out" 2> "$artifact/$count-$mode.err"
 done
done
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 /usr/bin/time -v -o "$artifact/132-structure-response14.time" \
 "$artifact/test_structure_tiles" --structure-response14-once 132 \
 > "$artifact/132-structure-response14.out" 2> "$artifact/132-structure-response14.err"
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 env XMVB_PHASE_TIMING=1 "$artifact/test_structure_tiles" --structure-response12 132 \
 > "$artifact/phases.out" 2> "$artifact/phases.err"
