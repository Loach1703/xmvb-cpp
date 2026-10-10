#!/usr/bin/env bash
#SBATCH --job-name=mnf2-storage
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=15-00:00:00
#SBATCH --output=mnf2-storage-%j.out
#SBATCH --error=mnf2-storage-%j.err
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
unset OMP_WAIT_POLICY GOMP_SPINCOUNT
root="$SLURM_SUBMIT_DIR/mnf2-storage-$SLURM_JOB_ID"
mkdir -p "$root"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
scontrol show job "$SLURM_JOB_ID" > "$root/allocation.txt"
cmake -S "$source_dir" -B "$source_dir/build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG -g1 -fno-omit-frame-pointer' \
    -DBUILD_TESTING=ON -DXMVB_BUILD_DEV_TARGETS=ON \
    -DCMAKE_PREFIX_PATH="$dependencies" \
    -DCMAKE_C_COMPILER="$dependencies/bin/cc" \
    -DCMAKE_CXX_COMPILER="$dependencies/bin/c++" \
    -DCMAKE_Fortran_COMPILER="$dependencies/bin/gfortran" \
    -DXMVB_ENABLE_EIGEN_BLAS_BACKEND=ON > "$root/build.log" 2>&1
cmake --build "$source_dir/build" --target xmvb test_active_two_electron_sparse \
    test_ri_ao_h1e_hvp test_orthogonal_direct_ci -j16 >> "$root/build.log" 2>&1
for threads in 1 8 32; do
    for test in test_active_two_electron_sparse test_ri_ao_h1e_hvp test_orthogonal_direct_ci; do
        srun --exact --cpus-per-task="$threads" --cpu-bind=cores \
            env OMP_NUM_THREADS="$threads" "$source_dir/build/src/$test" \
            > "$root/$test-$threads.out" 2>&1
    done
done
cp "$source_dir/build/src/xmvb-cpp.exe" "$root/current.exe"
cp "$baseline" "$root/baseline.exe"
variants=(baseline-8 current-8 baseline-16 current-16 baseline-32 current-32 baseline-passive-16 current-passive-16 baseline-1 current-1)
if [[ "${3:-}" == isolated ]]; then
    variants=(baseline-8 current-8 baseline-16 current-16 baseline-32 current-32 baseline-1 current-1)
elif [[ "${3:-}" == reverse ]]; then
    # Reverse matched pairs to check order/background-load effects, without
    # another costly single-thread calibration.
    variants=(current-32 baseline-32 current-16 baseline-16 current-8 baseline-8)
fi
for variant in "${variants[@]}"; do
    threads="${variant##*-}"
    exe="$root/current.exe"
    [[ "$variant" != baseline-* ]] || exe="$root/baseline.exe"
    launch=(env -u OMP_WAIT_POLICY -u GOMP_SPINCOUNT)
    [[ "$variant" != *-passive-* ]] || launch=(env OMP_WAIT_POLICY=PASSIVE GOMP_SPINCOUNT=0)
    mkdir -p "$root/$variant"
    cd "$root/$variant"
    cp "$source_dir/testdata/vbscf/MnF2.xmi" input.xmi
    sha256sum "$exe" input.xmi > hashes.txt
    srun --exact --cpus-per-task=32 --cpu-bind=verbose,cores \
        "${launch[@]}" OMP_NUM_THREADS="$threads" XMVB_PHASE_TIMING=1 \
        python3 "$source_dir/benchmarks/orbital-preconditioner/affinity_probe.py" \
        --threads "$threads" --compact --record affinity.json -- \
        /usr/bin/time -v -o resource.time "$exe" input.xmi \
        --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        > output.out 2> run.err
done
