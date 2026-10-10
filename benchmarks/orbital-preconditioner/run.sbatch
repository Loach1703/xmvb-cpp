#!/usr/bin/env bash
#SBATCH --job-name=orbital-projector-regression
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --array=0-7%45
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=run-%A_%a.out
#SBATCH --error=run-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
updated="$1"
baseline="$2"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 LC_ALL=C
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
index="$SLURM_ARRAY_TASK_ID"
cases=(241_vb FeCl2)
methods=(8 9)
variants=(baseline updated)
case_name="${cases[index / 4]}"
method="${methods[(index / 2) % 2]}"
variant="${variants[index % 2]}"
exe="$baseline/build/src/xmvb-cpp.exe"
if [[ "$variant" == updated ]]; then exe="$updated/build/src/xmvb-cpp.exe"; fi
run_dir="$SLURM_SUBMIT_DIR/run-$SLURM_ARRAY_JOB_ID/$case_name-iscf$method-$variant"
mkdir -p "$run_dir"
cd "$run_dir"
sed -E "s/iscf=[0-9]+/iscf=$method/" "$updated/testdata/vbscf/$case_name.xmi" > input.xmi
hostname > node.txt
lscpu > cpu.txt
sha256sum "$exe" input.xmi > hashes.txt
env | sort > environment.txt
status=0
srun --cpu-bind=cores /usr/bin/time -v -o resource.time \
  "$exe" input.xmi --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
  > output.out 2> run.err || status=$?
printf '%s\n' "$status" > exit_code.txt
exit "$status"
