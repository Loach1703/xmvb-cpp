#!/usr/bin/env bash
#SBATCH --job-name=gram-tiles
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=gram-tiles-%j.out
#SBATCH --error=gram-tiles-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
export OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
artifact="$SLURM_SUBMIT_DIR/gram-tiles-$SLURM_JOB_ID"
mkdir -p "$artifact"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
cmake --build "$SLURM_SUBMIT_DIR/build" --target test_structure_tiles -j8 > "$artifact/build.log" 2>&1
sha256sum "$SLURM_SUBMIT_DIR/build/src/test_structure_tiles" > "$artifact/binary.sha256"
sha256sum "$SLURM_SUBMIT_DIR/src/vbscf/structures/assembly/overlap_tiles.cpp" \
 "$SLURM_SUBMIT_DIR/src/vbscf/structures/orthogonal_ci/exterior_transform.cpp" \
 "$SLURM_SUBMIT_DIR/tests/vbscf/unit/test_structure_tiles.cpp" > "$artifact/sources.sha256"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
 bash -c 'grep Cpus_allowed_list /proc/self/status' > "$artifact/affinity.txt"
for count in 1 3 20 132; do
 srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o "$artifact/$count.time" \
  "$SLURM_SUBMIT_DIR/build/src/test_structure_tiles" --gram12 "$count" \
  > "$artifact/$count.out" 2> "$artifact/$count.err"
done
