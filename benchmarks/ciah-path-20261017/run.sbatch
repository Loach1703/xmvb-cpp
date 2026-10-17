#!/usr/bin/env bash
#SBATCH --job-name=ciah-path
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=0-2
#SBATCH --output=run-%A_%a.out
#SBATCH --error=run-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# -ge 1 && $# -le 2 && $SLURM_CPUS_PER_TASK == 32 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0
set -u
source_dir=$(realpath -e "$1")
preconditioner=${2:-hessian-diagonal}
case "$preconditioner" in hessian-diagonal|one-electron) ;; *) exit 64 ;; esac
baseline=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle
cases=(240 241 FeCl2 F2 MnF2)
case_name=${cases[$SLURM_ARRAY_TASK_ID]}
run="$SLURM_SUBMIT_DIR/run-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$run"
cp "$SLURM_SUBMIT_DIR/inputs/$case_name.xmi" "$run/input.xmi"
cp "$source_dir/source.sha256" "$source_dir/parent-revision.txt" "$run/"
ln -s "$source_dir/basis" "$run/basis"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="$baseline/lib:${LD_LIBRARY_PATH:-}"
export XMVB_CIAH_TRACE=scalar CIAH_PATH_EXECUTABLE="$source_dir/build/src/xmvb-cpp.exe"
unset XMVB_NEO_TRACE XMVB_NEO_KEYFRAME_AUDIT XMVB_NEO_ORBITAL_STEP_AUDIT
unset XMVB_CIAH_CANDIDATE_AUDIT XMVB_NEO_CANDIDATE_AUDIT XMVB_PHASE_TIMING
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_|XMVB_)' > environment.txt
sha256sum "$CIAH_PATH_EXECUTABLE" input.xmi "$baseline/lib/"* > hashes.txt
printf '%s\n' "$preconditioner" > preconditioner.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -c '
  grep Cpus_allowed_list /proc/self/status > affinity.txt
  exec /usr/bin/time -v -o resource.time "$CIAH_PATH_EXECUTABLE" input.xmi \
    --optimizer-backend ciah --eigensolver davidson --orbital-preconditioner "$1" \
    --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --ciah-trace ciah.tsv
' bash "$preconditioner" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
[[ $status == 0 ]] || exit "$status"
grep -q '^Termination reason *: ciah_projected_gradient_tolerance' output.out
