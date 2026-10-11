#!/usr/bin/env bash
#SBATCH --job-name=mnf2-clock-profile
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15-00:00:00
#SBATCH --output=mnf2-clock-%j.out
#SBATCH --error=mnf2-clock-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
profile_threads="${2:-16}"
if ! [[ "$profile_threads" =~ ^[0-9]+$ ]] ||
   ((profile_threads < 1 || profile_threads > SLURM_CPUS_PER_TASK)); then
    printf '%s\n' 'Profile thread count must fit the allocated cores.' >&2
    exit 64
fi
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS="$profile_threads" OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
unset OMP_WAIT_POLICY GOMP_SPINCOUNT
root="$SLURM_SUBMIT_DIR/mnf2-clock-$SLURM_JOB_ID"
mkdir -p "$root"
cd "$root"
hostname > node.txt
lscpu > cpu.txt
cp "$source_dir/build/src/xmvb-cpp.exe" current.exe
cp "$source_dir/testdata/vbscf/MnF2.xmi" input.xmi
sha256sum current.exe input.xmi > hashes.txt
# Probe both PMU and software-event permissions, without changing kernel policy.
srun --exact --cpus-per-task=1 --cpu-bind=cores \
    perf stat -e cycles:u,instructions:u -- true > perf-probe.out 2> perf-probe.err || true
srun --exact --cpus-per-task=1 --cpu-bind=cores \
    perf stat -e cpu-clock:u -- true > perf-software.out 2> perf-software.err || true
srun --exact --cpus-per-task=32 --cpu-bind=cores \
    python3 "$source_dir/benchmarks/orbital-preconditioner/affinity_probe.py" \
    --threads "$profile_threads" --compact --record affinity.json -- \
    /usr/bin/gprofng collect app -p on -o capture.er ./current.exe input.xmi \
    --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
    > output.out 2> run.err
for report in functions calltree threads; do
    srun --exact --cpus-per-task=1 --cpu-bind=cores \
        /usr/bin/gprofng display text -limit 40 "-$report" capture.er \
        > "$report.txt" 2> "$report.err"
done
srun --exact --cpus-per-task=1 --cpu-bind=cores \
    /usr/bin/gprofng display text -header -statistics capture.er > quality.txt 2> quality.err
