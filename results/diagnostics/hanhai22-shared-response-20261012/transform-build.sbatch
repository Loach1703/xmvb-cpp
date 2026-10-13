#!/usr/bin/env bash
#SBATCH --job-name=transform-build
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=transform-build-%j.out
#SBATCH --error=transform-build-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
artifact="$SLURM_SUBMIT_DIR/fecl2-ritz-cost-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
cp "$SLURM_SUBMIT_DIR/transform-source.patch" "$artifact/source.patch"
cp "$SLURM_SUBMIT_DIR/transform-source.txt" "$artifact/source.txt"
cmake --build "$source_dir/build" --target xmvb test_davidson test_keyframe test_ciah test_coupled_structure test_orthogonal_direct_ci test_sparse_orbital_quotient test_hessian_diagonal -j8 > "$artifact/build.log" 2>&1
for threads in 1 32; do
 OMP_NUM_THREADS=$threads srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores "$source_dir/build/src/test_orthogonal_direct_ci" > "$artifact/exterior-$threads.out" 2> "$artifact/exterior-$threads.err"
done
for test_name in test_davidson test_keyframe test_ciah test_coupled_structure test_sparse_orbital_quotient test_hessian_diagonal; do
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores "$source_dir/build/src/$test_name" > "$artifact/$test_name.out" 2> "$artifact/$test_name.err"
done
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/"
sha256sum "$artifact/xmvb-cpp.exe" > "$artifact/hash.txt"
grep -E 'CMAKE_BUILD_TYPE:|CMAKE_CXX_FLAGS(_RELEASE)?:' "$source_dir/build/CMakeCache.txt" > "$artifact/flags.txt"
