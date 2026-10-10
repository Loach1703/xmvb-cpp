#!/usr/bin/env bash
#SBATCH --job-name=neo-wait-policy
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --array=0-2%3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=15-00:00:00
#SBATCH --output=wait-policy-%A_%a.out
#SBATCH --error=wait-policy-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
cases=(241_vb 240_vb FeCl2)
case_name="${cases[SLURM_ARRAY_TASK_ID]}"
root="$SLURM_SUBMIT_DIR/wait-policy-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$root"
hostname > "$root/node.txt"
for policy in default passive; do
    run_dir="$root/$policy"
    mkdir -p "$run_dir"
    cd "$run_dir"
    cp "$source_dir/testdata/vbscf/$case_name.xmi" input.xmi
    unset OMP_WAIT_POLICY GOMP_SPINCOUNT
    if [[ "$policy" == passive ]]; then
        export OMP_WAIT_POLICY=PASSIVE GOMP_SPINCOUNT=0
    fi
    env | sort > environment.txt
    cat /proc/loadavg > load.txt
    exe="$source_dir/build/src/xmvb-cpp.exe"
    sha256sum "$exe" input.xmi > hashes.txt
    srun --exact --cpu-bind=cores env XMVB_PHASE_TIMING=1 \
        /usr/bin/time -v -o resource.time "$exe" input.xmi \
        --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        > output.out 2> run.err
done
