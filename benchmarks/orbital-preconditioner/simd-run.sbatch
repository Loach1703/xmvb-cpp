#!/usr/bin/env bash
#SBATCH --job-name=simd-shell-pair
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --array=0-4%5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15-00:00:00
#SBATCH --output=simd-run-%A_%a.out
#SBATCH --error=simd-run-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
baseline="$1"
updated="$2"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
cases=(241_vb 240_vb FeCl2 MnF2 csd_CERRAS_vb)
case_name="${cases[SLURM_ARRAY_TASK_ID]}"
root="$SLURM_SUBMIT_DIR/simd-run-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$root"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
env | sort > "$root/environment.txt"
for variant in baseline updated; do
    source_dir="$baseline"
    [[ "$variant" != updated ]] || source_dir="$updated"
    run_dir="$root/$variant"
    mkdir -p "$run_dir"
    cd "$run_dir"
    cp "$baseline/testdata/vbscf/$case_name.xmi" input.xmi
    exe="$source_dir/build/src/xmvb-cpp.exe"
    sha256sum "$exe" input.xmi > hashes.txt
    srun --exact --cpu-bind=cores env XMVB_PHASE_TIMING=1 \
        /usr/bin/time -v -o resource.time "$exe" input.xmi \
        --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        > output.out 2> run.err
done
