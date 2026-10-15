#!/usr/bin/env bash
#SBATCH --job-name=neo-residual
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=15-00:00:00
#SBATCH --output=neo-residual-%j.out
#SBATCH --error=neo-residual-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
binary="$2"
case_name="${3:-FeCl2}"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1 LC_ALL=C
export XMVB_NEO_TRACE=1 XMVB_PHASE_TIMING=1
unset OMP_WAIT_POLICY GOMP_SPINCOUNT XMVB_NEO_CANDIDATE_AUDIT
if [[ ${4:-} == audit ]]; then export XMVB_NEO_CANDIDATE_AUDIT=1; fi
root="$SLURM_SUBMIT_DIR/neo-residual-$SLURM_JOB_ID"
mkdir -p "$root"
cd "$root"
cp "/home/guqqgroup/taoxia/xmvb-runs/all-optimizers-c2776c4/neo/$case_name/input.xmi" input.xmi
sha256sum "$binary" input.xmi > hashes.txt
srun --exact --cpus-per-task=32 --cpu-bind=cores \
    python3 "$source_dir/benchmarks/orbital-preconditioner/affinity_probe.py" \
    --threads 32 --compact --record affinity.json -- \
    /usr/bin/time -v -o resource.time "$binary" input.xmi \
    --optimizer-backend neo --eigensolver davidson \
    --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
    --neo-trace trace.tsv > output.out 2> run.err
