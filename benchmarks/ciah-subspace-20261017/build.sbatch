#!/usr/bin/env bash
#SBATCH --job-name=ciah-subspace-build
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=build-%j.out
#SBATCH --error=build-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0 cmake/3.31.9
set -u
source_dir=$(realpath -e "$1")
deps=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
baseline=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle
build="$source_dir/build"
export LD_LIBRARY_PATH="$baseline/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
cmake -S "$source_dir" -B "$build" -DXMVB_BUILD_DEV_TARGETS=ON \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_CXX_FLAGS_RELEASE='-O3 -DNDEBUG -g1 -fno-omit-frame-pointer' \
  -DEIGEN3_ROOT_DIR="$deps" -DCINT_ROOT_DIR="$deps" -DLAPACK_ROOT_DIR="$deps" \
  -DCINT_LIB="$baseline/lib/libcint.so.6" \
  -DLAPACK_LIB="$baseline/lib/libopenblas.so.0" \
  -DLAPACKE_LIB="$baseline/lib/libopenblas.so.0" > configure.log 2>&1
cmake --build "$build" --target xmvb test_ciah test_ci_transport \
  test_selected_ci_transport test_keyframe test_coupled_structure \
  test_hessian_diagonal test_orbital_step_geometry -j16 > build.log 2>&1
ctest --test-dir "$build" --output-on-failure \
  -R '^(ciah|ci_transport|selected_ci_transport|keyframe|coupled_structure|hessian_diagonal|orbital_step_geometry)$' \
  > ctest.out 2>&1
sha256sum "$build/src/xmvb-cpp.exe" > binary.sha256
cp "$build/CMakeCache.txt" .
module list > modules.txt 2>&1
touch build.passed
