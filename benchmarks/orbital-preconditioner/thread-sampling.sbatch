#!/usr/bin/env bash
#SBATCH --job-name=241-thread-samples
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=15-00:00:00
#SBATCH --output=thread-sampling-%j.out
#SBATCH --error=thread-sampling-%j.err
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
root="$SLURM_SUBMIT_DIR/thread-sampling-$SLURM_JOB_ID"
mkdir -p "$root"
cd "$root"
hostname > node.txt
cp "$source_dir/testdata/vbscf/241_vb.xmi" input.xmi
exe="$source_dir/build/src/xmvb-cpp.exe"
sha256sum "$exe" input.xmi > hashes.txt
srun --exact --cpu-bind=cores /usr/bin/gprofng collect app -p on -o capture.er \
    "$exe" input.xmi --optimizer-backend neo \
    --gradient-tolerance 1e-3 --energy-tolerance 1e-7 > output.out 2> run.err
for report in functions calltree threads; do
    srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
        /usr/bin/gprofng display text -limit 40 "-$report" capture.er \
        > "$report.txt" 2> "$report.err"
done
srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
    /usr/bin/gprofng display text -header -statistics capture.er > quality.txt 2> quality.err
