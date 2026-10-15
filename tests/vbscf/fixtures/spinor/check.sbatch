#!/bin/bash
#SBATCH --job-name=2c-p0p2
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
snapshot=/home/guqqgroup/taoxia/work/xmvb-2c-p0p2-20261009
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export PATH="$dependencies/bin:$PATH" CONDA_PREFIX="$dependencies"
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export CCACHE_DIR=/home/guqqgroup/taoxia/.cache/ccache-xmvb-cpp
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$snapshot"
hostname
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON -DXMVB_BUILD_DEV_TARGETS=ON \
  -DCMAKE_PREFIX_PATH="$dependencies" \
  -DCMAKE_C_COMPILER="$dependencies/bin/cc" \
  -DCMAKE_CXX_COMPILER="$dependencies/bin/c++" \
  -DCMAKE_Fortran_COMPILER="$dependencies/bin/gfortran" \
  -DXMVB_ENABLE_EIGEN_BLAS_BACKEND=ON > configure.log 2>&1
cmake --build build --target test_spinor_strings test_spinor_projection \
  test_orthogonal_direct_ci -j16 > build.log 2>&1
ctest --test-dir build -R '^(spinor_strings|spinor_projection|orthogonal_direct_ci)$' \
  --output-on-failure > tests.log 2>&1
date -Is > finished.txt
