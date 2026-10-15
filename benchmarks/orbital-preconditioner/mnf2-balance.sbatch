#!/usr/bin/env bash
#SBATCH --job-name=mnf2-worker-balance
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=15-00:00:00
#SBATCH --output=mnf2-balance-%j.out
#SBATCH --error=mnf2-balance-%j.err
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
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
export LC_ALL=C OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
unset OMP_WAIT_POLICY GOMP_SPINCOUNT
root="$SLURM_SUBMIT_DIR/mnf2-balance-$SLURM_JOB_ID"
mkdir -p "$root"
cmake -S "$source_dir" -B "$source_dir/build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG -g1 -fno-omit-frame-pointer' \
    -DBUILD_TESTING=ON -DXMVB_BUILD_DEV_TARGETS=ON \
    -DCMAKE_PREFIX_PATH="$dependencies" \
    -DCMAKE_C_COMPILER="$dependencies/bin/cc" \
    -DCMAKE_CXX_COMPILER="$dependencies/bin/c++" \
    -DCMAKE_Fortran_COMPILER="$dependencies/bin/gfortran" \
    -DXMVB_ENABLE_EIGEN_BLAS_BACKEND=ON > "$root/build.log" 2>&1
cmake --build "$source_dir/build" --target xmvb test_active_two_electron_sparse \
    test_ri_ao_h1e_hvp -j16 >> "$root/build.log" 2>&1
for test in test_active_two_electron_sparse test_ri_ao_h1e_hvp; do
    srun --exact --cpus-per-task=32 --cpu-bind=cores \
        env OMP_NUM_THREADS=32 XMVB_PHASE_TIMING=1 "$source_dir/build/src/$test" \
        > "$root/$test.out" 2> "$root/$test.err"
done
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
cp "$source_dir/build/src/xmvb-cpp.exe" "$root/current.exe"
for threads in 8 16 32; do
    mkdir -p "$root/current-$threads"
    cd "$root/current-$threads"
    cp "$source_dir/testdata/vbscf/MnF2.xmi" input.xmi
    sha256sum "$root/current.exe" input.xmi > hashes.txt
    srun --exact --cpus-per-task=32 --cpu-bind=cores \
        env OMP_NUM_THREADS="$threads" XMVB_PHASE_TIMING=1 \
        python3 "$source_dir/benchmarks/orbital-preconditioner/affinity_probe.py" \
        --threads "$threads" --compact --record affinity.json -- \
        /usr/bin/time -v -o resource.time "$root/current.exe" input.xmi \
        --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        > output.out 2> run.err
done
