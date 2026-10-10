#!/usr/bin/env bash
#SBATCH --job-name=mnf2-metric
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15-00:00:00
#SBATCH --output=mnf2-metric-%j.out
#SBATCH --error=mnf2-metric-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
baseline="$2"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export PATH="$dependencies/bin:$PATH"
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export CCACHE_DIR=/home/guqqgroup/taoxia/.cache/ccache-xmvb-cpp
export LC_ALL=C OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
root="$SLURM_SUBMIT_DIR/mnf2-metric-$SLURM_JOB_ID"
mkdir -p "$root"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
cmake --build "$source_dir/build" --target xmvb test_sparse_orbital_quotient \
    test_coupled_structure test_shifted_metric_preconditioner -j16 > "$root/build.log" 2>&1
for threads in 1 8 32; do
    for test in test_sparse_orbital_quotient test_coupled_structure test_shifted_metric_preconditioner; do
        srun --exact --cpus-per-task="$threads" --cpu-bind=cores \
            env OMP_NUM_THREADS="$threads" "$source_dir/build/src/$test" \
            > "$root/$test-$threads.out" 2>&1
    done
done
for variant in baseline updated; do
    run_dir="$root/$variant"
    mkdir -p "$run_dir"
    cd "$run_dir"
    cp "$source_dir/testdata/vbscf/MnF2.xmi" input.xmi
    exe="$baseline"
    [[ "$variant" != updated ]] || exe="$source_dir/build/src/xmvb-cpp.exe"
    sha256sum "$exe" input.xmi > hashes.txt
    srun --exact --cpus-per-task=32 --cpu-bind=cores \
        env OMP_NUM_THREADS=32 XMVB_PHASE_TIMING=1 \
        /usr/bin/time -v -o resource.time "$exe" input.xmi \
        --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        > output.out 2> run.err
done
