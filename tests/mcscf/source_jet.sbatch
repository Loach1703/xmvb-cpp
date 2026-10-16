#!/usr/bin/env bash
#SBATCH --job-name=source-jet-check
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:20:00
#SBATCH --output=source-jet-%j.out
#SBATCH --error=source-jet-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
artifact="$source_dir/source-jet-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
module list > "$artifact/modules.txt" 2>&1
cmake --build "$source_dir/build" --target xmvb test_structure_tiles test_orthogonal_direct_ci \
 test_coupled_structure test_keyframe test_ciah test_sparse_orbital_quotient \
 test_support_preserving_gauge -j8 > "$artifact/build.log" 2>&1
sha256sum "$source_dir/build/src/xmvb-cpp.exe" "$source_dir/build/src/test_structure_tiles" > "$artifact/binaries.sha256"
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores \
 ctest --test-dir "$source_dir/build" --output-on-failure \
 -R '^(structure_tiles|orthogonal_direct_ci|coupled_structure|source_cover|keyframe|ciah|sparse_orbital_quotient|support_preserving_gauge)$' \
 > "$artifact/ctest.out" 2>&1
# CAS12 checks the scalar shortcut against an independent full-product oracle,
# and a nonclosed three-column subspace by directional finite differences.
for count in 1 3; do
 args=(--product-jet12 "$count")
 [[ $count == 1 ]] && args=(--source-state-reference 12)
 srun --exact --ntasks=1 --cpus-per-task=8 --cpu-bind=cores \
  python3 "$source_dir/tests/mcscf/pin_ci.py" 8 "$artifact/affinity.jsonl" \
  "$source_dir/build/src/test_structure_tiles" "${args[@]}" \
  > "$artifact/cas12-$count.out" 2> "$artifact/cas12-$count.err"
done
