#!/usr/bin/env bash
#SBATCH --job-name=neo-expansion-audit
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=15-00:00:00
#SBATCH --output=neo-expansion-%j.out
#SBATCH --error=neo-expansion-%j.err
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
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1 LC_ALL=C
unset OMP_WAIT_POLICY GOMP_SPINCOUNT XMVB_NEO_TRACE
root="$SLURM_SUBMIT_DIR/neo-expansion-$SLURM_JOB_ID"
mkdir -p "$root"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
cmake --build "$source_dir/build" --target xmvb test_neo \
    test_neo_response_solver test_neo_globalization test_sparse_orbital_quotient \
    test_active_two_electron_sparse test_ri_ao_h1e_hvp \
    test_shifted_metric_preconditioner -j16 \
    > "$root/build.log" 2>&1
for test in test_neo test_neo_response_solver test_neo_globalization \
    test_sparse_orbital_quotient test_shifted_metric_preconditioner \
    test_active_two_electron_sparse test_ri_ao_h1e_hvp; do
    srun --exact --cpus-per-task=32 --cpu-bind=cores \
        "$source_dir/build/src/$test" > "$root/$test.out" 2> "$root/$test.err"
done
cp "$source_dir/build/src/xmvb-cpp.exe" "$root/current.exe"
for mode in native trace; do
    mkdir -p "$root/$mode"
    cd "$root/$mode"
    cp "$source_dir/testdata/vbscf/MnF2.xmi" input.xmi
    sha256sum "$root/current.exe" input.xmi > hashes.txt
    if [[ "$mode" == trace ]]; then export XMVB_NEO_TRACE=1; fi
    srun --exact --cpus-per-task=32 --cpu-bind=cores \
        python3 "$source_dir/benchmarks/orbital-preconditioner/affinity_probe.py" \
        --threads 32 --compact --record affinity.json -- \
        /usr/bin/time -v -o resource.time "$root/current.exe" input.xmi \
        --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        --neo-trace trace.tsv > output.out 2> run.err
done
