#!/usr/bin/env bash
#SBATCH --job-name=casscf-build
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
[[ -n ${SLURM_JOB_ID:-} && $# -ge 1 && $# -le 2 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0 cmake/3.31.9
set -u
source_dir=$(realpath -e "$1")
build_dir=${2:-$source_dir/build}
deps=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
bundle=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle
export LD_LIBRARY_PATH="$bundle/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
cmake -S "$source_dir" -B "$build_dir" -DXMVB_BUILD_DEV_TARGETS=ON \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_CXX_FLAGS_RELEASE='-O3 -DNDEBUG -g1 -fno-omit-frame-pointer' \
  -DEIGEN3_ROOT_DIR="$deps" -DCINT_ROOT_DIR="$deps" -DLAPACK_ROOT_DIR="$deps" \
  -DCINT_LIB="$bundle/lib/libcint.so.6" \
  -DLAPACK_LIB="$bundle/lib/libopenblas.so.0" \
  -DLAPACKE_LIB="$bundle/lib/libopenblas.so.0" > configure.log 2>&1
cmake --build "$build_dir" --target xmvb casscf_fixture test_casci \
  test_casscf_derivatives test_casscf_coupled test_ciah test_ci_transport \
  test_selected_ci_transport test_keyframe test_coupled_structure \
  test_hessian_diagonal test_orbital_step_geometry test_input_guess_selection -j16 > build.log 2>&1
ctest --test-dir "$build_dir" --output-on-failure \
  -R '^(casci|casscf_xmi|casscf_derivatives|casscf_coupled|ciah|ci_transport|selected_ci_transport|keyframe|coupled_structure|hessian_diagonal|orbital_step_geometry|input_guess_selection)$' \
  > ctest.out 2>&1
sha256sum "$build_dir/src/"{xmvb-cpp.exe,casscf_fixture} > binary.sha256
cp "$build_dir/CMakeCache.txt" .
module list > modules.txt 2>&1
