#!/usr/bin/env bash
#SBATCH --job-name=neo-history-panel
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --array=0-21
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=15-00:00:00
#SBATCH --output=neo-panel-%A_%a.out
#SBATCH --error=neo-panel-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
snapshot="$2"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1 LC_ALL=C
unset OMP_WAIT_POLICY GOMP_SPINCOUNT XMVB_NEO_TRACE
mapfile -t cases < "$source_dir/benchmarks/orbital-preconditioner/neo-panel-cases.txt"
case_name="${cases[SLURM_ARRAY_TASK_ID]}"
historical="/home/guqqgroup/taoxia/xmvb-runs/all-optimizers-c2776c4/neo/$case_name"
root="$SLURM_SUBMIT_DIR/neo-panel-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$root/historical" "$root/current"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
module list > "$root/modules.txt" 2>&1
for variable in OMP_NUM_THREADS OMP_PLACES OMP_PROC_BIND OMP_MAX_ACTIVE_LEVELS \
    OPENBLAS_NUM_THREADS MKL_NUM_THREADS BLIS_NUM_THREADS; do
    printf '%s=%s\n' "$variable" "${!variable}" >> "$root/threads.txt"
done
cp "$historical/input.xmi" "$historical/optimizer.tsv" \
    "$historical/run.out" "$historical/resource.time" "$root/historical/"
sha256sum "$historical/input.xmi" > "$root/historical/input.sha256"
cp "$snapshot" "$root/current.exe"
cd "$root/current"
cp "$historical/input.xmi" input.xmi
sha256sum "$root/current.exe" input.xmi > hashes.txt
set +e
srun --exact --cpus-per-task=32 --cpu-bind=cores \
    python3 "$source_dir/benchmarks/orbital-preconditioner/affinity_probe.py" \
    --threads 32 --compact --record affinity.json -- \
    /usr/bin/time -v -o resource.time "$root/current.exe" input.xmi \
    --optimizer-backend neo --eigensolver davidson \
    --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
    --neo-trace trace.tsv > output.out 2> run.err
run_status=$?
printf '%s\n' "$run_status" > status.txt
exit "$run_status"
