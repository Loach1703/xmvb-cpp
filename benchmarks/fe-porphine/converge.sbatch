#!/usr/bin/env bash
#SBATCH --job-name=fe-corrected-guess
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --array=0-5%6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=15-00:00:00
#SBATCH --output=converge-%A_%a.out
#SBATCH --error=converge-%A_%a.err

set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
prepared="$1"
source_dir="$2"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export PATH="$dependencies/bin:$PATH" CONDA_PREFIX="$dependencies"
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
states=(singlet triplet quintet)
methods=(8 9)
index="$SLURM_ARRAY_TASK_ID"
[[ "$index" =~ ^[0-5]$ ]] || exit 64
state="${states[index / 2]}"
iscf="${methods[index % 2]}"
exe="$source_dir/build/src/xmvb-cpp.exe"
[[ -x "$exe" && -x /usr/bin/time && -f "$prepared/$state-iscf$iscf.xmi" ]] || exit 66
run_dir="$SLURM_SUBMIT_DIR/converge-$SLURM_ARRAY_JOB_ID/$state-iscf$iscf"
mkdir -p "$run_dir"
cp "$prepared/$state-iscf$iscf.xmi" "$run_dir/input.xmi"
cp "$prepared/initial_audit.json" "$prepared/casci-check.json" "$run_dir/"
cp "$SLURM_SUBMIT_DIR/source-revision.txt" "$SLURM_SUBMIT_DIR/converge.sbatch" "$run_dir/"
cd "$run_dir"
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
for variable in OMP_NUM_THREADS OMP_PLACES OMP_PROC_BIND OMP_MAX_ACTIVE_LEVELS \
    OPENBLAS_NUM_THREADS MKL_NUM_THREADS BLIS_NUM_THREADS GOTO_NUM_THREADS \
    SLURM_JOB_ID SLURM_ARRAY_JOB_ID SLURM_ARRAY_TASK_ID SLURM_CPUS_PER_TASK; do
    printf '%s=%s\n' "$variable" "${!variable-}"
done > allocation.txt
sha256sum "$exe" input.xmi "$prepared/initial_mo.npy" > hashes.txt
date -Is > started.txt
status=0
# No shell timeout. Only the mandatory partition maximum (15 days) applies.
# Retain the input's 1000-iteration safety limit, not a shortened benchmark cap.
srun --cpu-bind=cores /usr/bin/time -v -o resource.time \
    "$exe" input.xmi --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
    > output.out 2> run.err || status=$?
printf 'state,iscf,exit_code\n%s,%s,%s\n' "$state" "$iscf" "$status" > exit_codes.csv
date -Is > finished.txt
exit "$status"
