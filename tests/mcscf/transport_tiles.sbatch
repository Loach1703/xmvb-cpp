#!/usr/bin/env bash
#SBATCH --job-name=transport-tiles
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:15:00
#SBATCH --output=transport-tiles-%j.out
#SBATCH --error=transport-tiles-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
artifact="$source_dir/transport-tiles-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
module list > "$artifact/modules.txt" 2>&1
cmake --build "$source_dir/build" --target xmvb test_representation_overlap \
 test_ci_transport test_selected_ci_transport test_keyframe test_structure_tiles \
 test_sparse_orbital_quotient -j8 > "$artifact/build.log" 2>&1
sha256sum "$source_dir/build/src/xmvb-cpp.exe" \
 "$source_dir/build/src/test_representation_overlap" > "$artifact/binaries.sha256"
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
srun --exact --ntasks=1 --cpus-per-task=4 --cpu-bind=cores \
 ctest --test-dir "$source_dir/build" --output-on-failure \
 -R '^(representation_overlap|ci_transport|selected_ci_transport|structure_tiles|keyframe|sparse_orbital_quotient|neo_hessian_diagonal_f2_hao|neo_hessian_diagonal_f2_oeo)$' \
 > "$artifact/ctest.out" 2>&1
cp "$source_dir/build/Testing/Temporary/LastTest.log" "$artifact/test-details.log"
