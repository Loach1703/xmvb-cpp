#!/usr/bin/env bash
#SBATCH --job-name=neo-admission-panel
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --array=0-4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=15-00:00:00
#SBATCH --output=neo-admission-panel-%A_%a.out
#SBATCH --error=neo-admission-panel-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
snapshots="$2"
candidate="${3:-coupled}"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1 LC_ALL=C
unset OMP_WAIT_POLICY GOMP_SPINCOUNT XMVB_NEO_TRACE
cases=(240_vb 67970_vb csd_CONYAH_vb csd_QOCGUM_vb csd_YEPKAF_vb)
case_name="${cases[SLURM_ARRAY_TASK_ID]}"
historical="/home/guqqgroup/taoxia/xmvb-runs/all-optimizers-c2776c4/neo/$case_name"
root="$SLURM_SUBMIT_DIR/neo-admission-panel-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$root"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
module list > "$root/modules.txt" 2>&1
for variable in OMP_NUM_THREADS OMP_PLACES OMP_PROC_BIND OMP_MAX_ACTIVE_LEVELS \
    OPENBLAS_NUM_THREADS MKL_NUM_THREADS BLIS_NUM_THREADS; do
    printf '%s=%s\n' "$variable" "${!variable}" >> "$root/threads.txt"
done
# Independent systems run in parallel; each matched pair shares one node.
for variant in current "$candidate"; do
    mkdir -p "$root/$variant"
    cd "$root/$variant"
    cp "$historical/input.xmi" input.xmi
    sha256sum "$snapshots/$variant.exe" input.xmi > hashes.txt
    set +e
    srun --exact --cpus-per-task=32 --cpu-bind=cores \
        python3 "$source_dir/benchmarks/orbital-preconditioner/affinity_probe.py" \
        --threads 32 --compact --record affinity.json -- \
        /usr/bin/time -v -o resource.time "$snapshots/$variant.exe" input.xmi \
        --optimizer-backend neo --eigensolver davidson \
        --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        --neo-trace trace.tsv > output.out 2> run.err
    run_status=$?
    set -e
    printf '%s\n' "$run_status" > status.txt
done
