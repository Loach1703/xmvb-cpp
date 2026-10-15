#!/usr/bin/env bash
#SBATCH --job-name=native-openmp-check
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=15-00:00:00
#SBATCH --output=native-check-%j.out
#SBATCH --error=native-check-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
baseline="$2"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export PATH="$dependencies/bin:$PATH" CONDA_PREFIX="$dependencies"
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export CCACHE_DIR=/home/guqqgroup/taoxia/.cache/ccache-xmvb-cpp
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
export LC_ALL=C OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
root="$SLURM_SUBMIT_DIR/native-check-$SLURM_JOB_ID"
mkdir -p "$root"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
module list > "$root/modules.txt" 2>&1
cmake -S "$source_dir" -B "$source_dir/build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG -g1 -fno-omit-frame-pointer' \
    -DBUILD_TESTING=ON -DXMVB_BUILD_DEV_TARGETS=ON \
    -DCMAKE_PREFIX_PATH="$dependencies" \
    -DCMAKE_C_COMPILER="$dependencies/bin/cc" \
    -DCMAKE_CXX_COMPILER="$dependencies/bin/c++" \
    -DCMAKE_Fortran_COMPILER="$dependencies/bin/gfortran" \
    -DXMVB_ENABLE_EIGEN_BLAS_BACKEND=ON > "$root/build.log" 2>&1
cmake --build "$source_dir/build" --target xmvb test_libcint_ri \
    test_active_two_electron_sparse test_orthogonal_direct_ci -j16 >> "$root/build.log" 2>&1
for threads in 1 8 32; do
    for test in test_libcint_ri test_active_two_electron_sparse test_orthogonal_direct_ci; do
        srun --exact --cpus-per-task="$threads" --cpu-bind=cores \
            env OMP_NUM_THREADS="$threads" "$source_dir/build/src/$test" \
            > "$root/$test-$threads.out" 2>&1
    done
done
cp "$source_dir/build/src/xmvb-cpp.exe" "$root/current.exe"
cp "$baseline" "$root/baseline.exe"
sha256sum "$root"/*.exe > "$root/binaries.sha256"
for variant in baseline-32 current-8 current-16 current-32; do
    threads="${variant##*-}"
    exe="$root/current.exe"
    [[ "$variant" != baseline-32 ]] || exe="$root/baseline.exe"
    mkdir -p "$root/$variant"
    cd "$root/$variant"
    cp "$source_dir/testdata/vbscf/MnF2.xmi" input.xmi
    sha256sum input.xmi > input.sha256
    srun --exact --cpus-per-task="$threads" --cpu-bind=cores \
        env OMP_NUM_THREADS="$threads" XMVB_PHASE_TIMING=1 \
        /usr/bin/time -v -o resource.time "$exe" input.xmi \
        --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        > output.out 2> run.err
done
