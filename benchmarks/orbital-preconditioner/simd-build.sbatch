#!/usr/bin/env bash
#SBATCH --job-name=simd-shell-build
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=15-00:00:00
#SBATCH --output=simd-build-%j.out
#SBATCH --error=simd-build-%j.err
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
root="$SLURM_SUBMIT_DIR/simd-build-$SLURM_JOB_ID"
mkdir -p "$root"
cd "$root"
cmake -S "$source_dir" -B "$source_dir/build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG -g1 -fno-omit-frame-pointer' \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    -DBUILD_TESTING=ON -DXMVB_BUILD_DEV_TARGETS=ON \
    -DCMAKE_PREFIX_PATH="$dependencies" \
    -DCMAKE_C_COMPILER="$dependencies/bin/cc" \
    -DCMAKE_CXX_COMPILER="$dependencies/bin/c++" \
    -DCMAKE_Fortran_COMPILER="$dependencies/bin/gfortran" \
    -DXMVB_ENABLE_EIGEN_BLAS_BACKEND=ON > build.log 2>&1
cmake --build "$source_dir/build" --target xmvb test_libcint_ri \
    test_orthogonal_direct_ci check_direct_libcint_smoke check_libcint_ri_smoke \
    -j16 >> build.log 2>&1
for threads in 1 8 32; do
    for test in test_libcint_ri test_orthogonal_direct_ci; do
        srun --exact --cpus-per-task="$threads" --cpu-bind=cores \
            env OMP_NUM_THREADS="$threads" "$source_dir/build/src/$test" \
            > "$test-$threads.out" 2>&1
    done
done
python3 "$source_dir/benchmarks/orbital-preconditioner/vector-report.py" \
    "$source_dir/build/compile_commands.json" "$root"
