#!/usr/bin/env bash
#SBATCH --job-name=neo-admission
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=15-00:00:00
#SBATCH --output=neo-admission-%j.out
#SBATCH --error=neo-admission-%j.err
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
root="$SLURM_SUBMIT_DIR/neo-admission-$SLURM_JOB_ID"
mkdir -p "$root"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
backend="$source_dir/src/vbscf/optimization/backends/neo.cpp"
executable="$source_dir/build/src/xmvb-cpp.exe"
cp "$backend" "$root/original-neo.cpp"
cp "$executable" "$root/current.exe"
sha256sum "$backend" "$executable" > "$root/original.sha256"
cp "$source_dir/benchmarks/orbital-preconditioner/neo-admission.patch" "$root/"
# Only the backend selection changes; restore mirror source and executable
# even if compilation or a run fails. Touch source to invalidate patched objects.
restore_mirror() {
    cp "$root/original-neo.cpp" "$backend"
    touch "$backend"
    cp "$root/current.exe" "$executable"
}
trap restore_mirror EXIT
git -C "$source_dir" apply --check "$root/neo-admission.patch"
git -C "$source_dir" apply "$root/neo-admission.patch"
cmake --build "$source_dir/build" --target xmvb -j16 > "$root/build.log" 2>&1
cp "$executable" "$root/coupled.exe"
restore_mirror
trap - EXIT
historical=/home/guqqgroup/taoxia/xmvb-runs/all-optimizers-c2776c4/neo/csd_NEGPIA_vb
# Sequential variants intentionally share the same allocation and affinity.
# They are a controlled comparison, not independent dataset jobs.
for variant in current coupled; do
    mkdir -p "$root/$variant"
    cd "$root/$variant"
    cp "$historical/input.xmi" input.xmi
    sha256sum "$root/$variant.exe" input.xmi > hashes.txt
    set +e
    srun --exact --cpus-per-task=32 --cpu-bind=cores \
        python3 "$source_dir/benchmarks/orbital-preconditioner/affinity_probe.py" \
        --threads 32 --compact --record affinity.json -- \
        /usr/bin/time -v -o resource.time "$root/$variant.exe" input.xmi \
        --optimizer-backend neo --eigensolver davidson \
        --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        --neo-trace trace.tsv > output.out 2> run.err
    run_status=$?
    set -e
    printf '%s\n' "$run_status" > status.txt
done
