#!/usr/bin/env bash
#SBATCH --job-name=master-regression
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=02:00:00
#SBATCH --output=master-regression-%A_%a.out
#SBATCH --error=master-regression-%A_%a.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
cases=(F2-OEO FeCl2 MnF2 241)
case_name="${cases[$SLURM_ARRAY_TASK_ID]}"
run="$root/master-regression-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$run"
case "$case_name" in
  F2-OEO|FeCl2|MnF2) fixture="$root/ci-ray-regression-1252666/$case_name" ;;
  *) fixture="$root/inputs/$case_name" ;;
esac
cp "$fixture/input.xmi" "$fixture/"*.gbs "$run/"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export XMVB_NEO_TRACE=scalar
unset XMVB_NEO_KEYFRAME_AUDIT XMVB_NEO_ORBITAL_STEP_AUDIT XMVB_NEO_ORBITAL_PATH
binary="/home/guqqgroup/xiatao/xmvb-runs/pyrrole-recheck-20261014/horizontal-build-1255480/xmvb-cpp.exe"
printf '%s  %s\n' a5313100b12fd57e7d876e0d21c9333df5cc9d1a4d4dfca126c25b3f5da9b602 "$binary" | sha256sum -c -
arguments=(--optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal
  --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 --neo-trace neo.tsv)
hostname > node.txt
lscpu > cpu.txt
printf '%s\n' 'master source 39201b81; build 1255480; uncommitted RI edits excluded' > source.txt
env | sort > environment.txt
sha256sum "$binary" input.xmi *.gbs > hashes.txt
printf '%s\n' "${arguments[@]}" > arguments.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time "$binary" input.xmi "${arguments[@]}" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
