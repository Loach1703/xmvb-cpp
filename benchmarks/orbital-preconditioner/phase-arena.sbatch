#!/usr/bin/env bash
#SBATCH --job-name=orbital-arena-check
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --array=0-1%45
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --mem=96G
#SBATCH --time=00:10:00
#SBATCH --output=phase-arena-%A_%a.out
#SBATCH --error=phase-arena-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_root="$1"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export PATH="/home/guqqgroup/taoxia/software/bin:$dependencies/bin:$PATH"
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 LC_ALL=C
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
cases=(241_vb FeCl2)
methods=(8 9)
case_name="${cases[SLURM_ARRAY_TASK_ID]}"
method="${methods[SLURM_ARRAY_TASK_ID]}"
root="$SLURM_SUBMIT_DIR/phase-arena-$SLURM_ARRAY_JOB_ID/$case_name-iscf$method"
mkdir -p "$root"
hostname > "$root/node.txt"
env | sort > "$root/environment.txt"
srun --exact --cpu-bind=cores heaptrack_print --help > "$root/heaptrack-print-help.txt" 2>&1
# Diagnostic allocator ablation only; do not change application defaults.
for variant in baseline updated; do
  run_dir="$root/$variant"
  mkdir -p "$run_dir"
  cd "$run_dir"
  sed -E "s/iscf=[0-9]+/iscf=$method/" "$source_root/updated/testdata/vbscf/$case_name.xmi" > input.xmi
  exe="$source_root/$variant/build/src/xmvb-cpp.exe"
  sha256sum "$exe" input.xmi > hashes.txt
  srun --exact --cpu-bind=cores env XMVB_PHASE_TIMING=1 MALLOC_ARENA_MAX=1 \
    /usr/bin/time -v -o resource.time "$exe" input.xmi \
    --gradient-tolerance 1e-3 --energy-tolerance 1e-7 > output.out 2> run.err
done
