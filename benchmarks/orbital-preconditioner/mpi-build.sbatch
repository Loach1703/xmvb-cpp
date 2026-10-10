#!/usr/bin/env bash
#SBATCH --job-name=mpi-integral-build
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=15-00:00:00
#SBATCH --output=mpi-build-%j.out
#SBATCH --error=mpi-build-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06 openmpi/5.0.9
set -u
source_dir="$1"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export PATH="$dependencies/bin:$PATH" CONDA_PREFIX="$dependencies"
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export CCACHE_DIR=/home/guqqgroup/taoxia/.cache/ccache-xmvb-cpp
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
cmake -S "$source_dir" -B "$source_dir/build" -G Ninja \
    '-UMPI_*' \
    -DCMAKE_BUILD_TYPE=Release '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG -g1 -fno-omit-frame-pointer' \
    -DBUILD_TESTING=ON -DXMVB_BUILD_DEV_TARGETS=ON -DXMVB_BUILD_MPI_CHECK=ON \
    -DCMAKE_PREFIX_PATH="$dependencies" \
    -DCMAKE_C_COMPILER="$dependencies/bin/cc" \
    -DCMAKE_CXX_COMPILER="$dependencies/bin/c++" \
    -DCMAKE_Fortran_COMPILER="$dependencies/bin/gfortran" \
    -DMPI_CXX_COMPILER=/opt/openmpi/5.0.9-ucx/bin/mpicxx \
    '-DCMAKE_EXE_LINKER_FLAGS=-Wl,--allow-shlib-undefined -Wl,-rpath-link,/usr/lib/x86_64-linux-gnu -Wl,-rpath-link,/usr/lpp/mmfs/lib'
cmake --build "$source_dir/build" --target xmvb check_integral_shards \
    test_libcint_ri test_active_two_electron_sparse -j16
for threads in 1 8 32; do
    for test in test_libcint_ri test_active_two_electron_sparse; do
        srun --exact --cpus-per-task="$threads" --cpu-bind=cores \
            env OMP_NUM_THREADS="$threads" "$source_dir/build/src/$test"
    done
done
