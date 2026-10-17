#!/usr/bin/env bash
#SBATCH --job-name=ciah-subspace-regression
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --array=0-3
#SBATCH --output=regression-%A_%a.out
#SBATCH --error=regression-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0
set -u
source_dir=$(realpath -e "$1")
baseline=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle
cases=(F2 241_vb FeCl2 241_reference)
case_name=${cases[$SLURM_ARRAY_TASK_ID]}
run="$SLURM_SUBMIT_DIR/regression-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$run"
if [[ $case_name == 241_reference ]]; then
  cp "$baseline/smoke/input.xmi" "$run/input.xmi"
else
  cp "$source_dir/testdata/vbscf/$case_name.xmi" "$run/input.xmi"
fi
ln -s "$source_dir/basis" "$run/basis"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="$baseline/lib:${LD_LIBRARY_PATH:-}"
export XMVB_CIAH_TRACE=scalar CIAH_SUBSPACE_EXECUTABLE="$source_dir/build/src/xmvb-cpp.exe"
hostname > node.txt
lscpu > cpu.txt
sha256sum "$CIAH_SUBSPACE_EXECUTABLE" input.xmi > hashes.txt
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -c '
  grep Cpus_allowed_list /proc/self/status > affinity.txt
  exec /usr/bin/time -v -o resource.time "$CIAH_SUBSPACE_EXECUTABLE" input.xmi \
    --optimizer-backend ciah --eigensolver davidson --orbital-preconditioner hessian-diagonal \
    --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --ciah-trace ciah.tsv
' > output.out 2> run.err
grep -q '^Termination reason *: ciah_projected_gradient_tolerance' output.out
