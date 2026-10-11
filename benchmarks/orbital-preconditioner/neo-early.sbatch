#!/usr/bin/env bash
#SBATCH --job-name=neo-early-pair
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --array=0-4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15-00:00:00
#SBATCH --output=neo-early-%A_%a.out
#SBATCH --error=neo-early-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
baseline="$2"
gradient_tolerance="${3:-1e-3}"
energy_tolerance="${4:-1e-7}"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1 LC_ALL=C
unset OMP_WAIT_POLICY GOMP_SPINCOUNT XMVB_NEO_TRACE
cases=(241_vb 240_vb FeCl2 MnF2 F2)
case_name="${cases[SLURM_ARRAY_TASK_ID]}"
root="$SLURM_SUBMIT_DIR/neo-early-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$root"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
cp "$source_dir/build/src/xmvb-cpp.exe" "$root/updated.exe"
for variant in baseline updated; do
    exe="$baseline"
    [[ "$variant" != updated ]] || exe="$root/updated.exe"
    mkdir -p "$root/$variant"
    cd "$root/$variant"
    cp "$source_dir/testdata/vbscf/$case_name.xmi" input.xmi
    sha256sum "$exe" input.xmi > hashes.txt
    srun --exact --cpus-per-task=32 --cpu-bind=cores \
        python3 "$source_dir/benchmarks/orbital-preconditioner/affinity_probe.py" \
        --threads 32 --compact --record affinity.json -- \
        /usr/bin/time -v -o resource.time "$exe" input.xmi \
        --optimizer-backend neo --gradient-tolerance "$gradient_tolerance" \
        --energy-tolerance "$energy_tolerance" \
        --neo-trace trace.tsv > output.out 2> run.err
done
