#!/usr/bin/env bash
#SBATCH --job-name=neo-diagonal-check
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=15-00:00:00
#SBATCH --output=neo-diagonal-check-%j.out
#SBATCH --error=neo-diagonal-check-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
binary="$2"
reference="${3:-}"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1 LC_ALL=C
unset XMVB_NEO_TRACE XMVB_PHASE_TIMING OMP_WAIT_POLICY GOMP_SPINCOUNT
root="$SLURM_SUBMIT_DIR/neo-diagonal-check-$SLURM_JOB_ID"
mkdir -p "$root"
failed=0
for case_name in F2 F2_OEO F2_SA2; do
    for mode in exact ri; do
      for variant in candidate reference; do
        run_binary="$binary"
        if [[ $variant == reference ]]; then
            [[ -n $reference ]] || continue
            run_binary="$reference"
        fi
        mkdir -p "$root/$case_name-$mode/$variant"
        cd "$root/$case_name-$mode/$variant"
        cp "$source_dir/testdata/vbscf/$case_name.xmi" input.xmi
        sha256sum "$run_binary" input.xmi > hashes.txt
        set +e
        srun --exact --cpus-per-task=32 --cpu-bind=cores \
            python3 "$source_dir/benchmarks/orbital-preconditioner/affinity_probe.py" \
            --threads 32 --compact --record affinity.json -- \
            "$run_binary" input.xmi --optimizer-backend neo --eigensolver davidson \
            --standard-two-electron-mode "$mode" \
            --orbital-preconditioner hessian-diagonal \
            --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
            --neo-trace trace.tsv > output.out 2> run.err
        run_status=$?
        set -e
        printf '%s\n' "$run_status" > status.txt
        if (( run_status != 0 )); then failed=2; fi
      done
    done
done
exit "$failed"
