#!/usr/bin/env bash
#SBATCH --job-name=neo-quality
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15-00:00:00
#SBATCH --output=neo-quality-%A_%a.out
#SBATCH --error=neo-quality-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1 LC_ALL=C
unset OMP_WAIT_POLICY GOMP_SPINCOUNT XMVB_NEO_TRACE XMVB_NEO_CANDIDATE_AUDIT
cases=(MnF2 FeCl2)
case_name="${cases[SLURM_ARRAY_TASK_ID]}"
root="$SLURM_SUBMIT_DIR/neo-quality-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$root"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
for mode in native audit; do
    mkdir -p "$root/$mode"
    cd "$root/$mode"
    cp "$source_dir/testdata/vbscf/$case_name.xmi" input.xmi
    sha256sum "$source_dir/build/src/xmvb-cpp.exe" input.xmi > hashes.txt
    if [[ "$mode" == audit ]]; then
        export XMVB_NEO_TRACE=1 XMVB_NEO_CANDIDATE_AUDIT=1
    fi
    srun --exact --cpus-per-task=32 --cpu-bind=cores \
        python3 "$source_dir/benchmarks/orbital-preconditioner/affinity_probe.py" \
        --threads 32 --compact --record affinity.json -- \
        /usr/bin/time -v -o resource.time "$source_dir/build/src/xmvb-cpp.exe" input.xmi \
        --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        --neo-trace trace.tsv > output.out 2> run.err
done
