#!/usr/bin/env bash
#SBATCH --job-name=orbital-projector-build
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=build-%j.out
#SBATCH --error=build-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export PATH="$dependencies/bin:$PATH" CONDA_PREFIX="$dependencies"
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export CCACHE_DIR=/home/guqqgroup/taoxia/.cache/ccache-xmvb-cpp
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1 LC_ALL=C
run_dir="$SLURM_SUBMIT_DIR/build-$SLURM_JOB_ID"
mkdir -p "$run_dir"
cd "$run_dir"
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
cmake -S "$source_dir" -B "$source_dir/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DXMVB_BUILD_DEV_TARGETS=ON \
  -DCMAKE_PREFIX_PATH="$dependencies" \
  -DCMAKE_C_COMPILER="$dependencies/bin/cc" \
  -DCMAKE_CXX_COMPILER="$dependencies/bin/c++" \
  -DCMAKE_Fortran_COMPILER="$dependencies/bin/gfortran" \
  -DXMVB_ENABLE_EIGEN_BLAS_BACKEND=ON > build.log 2>&1
cmake --build "$source_dir/build" --target xmvb test_projected_orbital_surrogate \
  test_normalized_orbital_curvature test_sparse_orbital_quotient test_neo_globalization \
  -j16 >> build.log 2>&1
for test in test_projected_orbital_surrogate test_normalized_orbital_curvature \
            test_sparse_orbital_quotient test_neo_globalization; do
  srun --cpu-bind=cores timeout 120s "$source_dir/build/src/$test" > "$test.out" 2>&1
done
srun --cpu-bind=cores /usr/bin/time -v -o kernel.time \
  "$source_dir/build/src/test_projected_orbital_surrogate" --benchmark > kernel.out 2>&1
sha256sum "$source_dir/build/src/xmvb-cpp.exe" > binary.sha256
