#!/usr/bin/env bash
#SBATCH --job-name=orbital-step-geometry
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=02:00:00
#SBATCH --output=step-geometry-%A_%a.out
#SBATCH --error=step-geometry-%A_%a.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
cases=(FeCl2 singlet triplet)
case_name="${cases[$SLURM_ARRAY_TASK_ID]}"
run="$root/step-geometry-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$run"
if [[ "$case_name" == FeCl2 ]]; then
  fixture="$root/ci-ray-regression-1252666/FeCl2"
else
  fixture="$root/inputs/$case_name"
  cap=2
fi
cp "$fixture/input.xmi" "$fixture/"*.gbs "$run/"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export XMVB_NEO_TRACE=scalar XMVB_NEO_ORBITAL_STEP_AUDIT=1
unset XMVB_NEO_KEYFRAME_AUDIT XMVB_NEO_CANDIDATE_AUDIT
binary="$root/step-geometry-build-${BUILD_JOB:?}/xmvb-cpp.exe"
arguments=(--optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal
  --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256
  --neo-trace neo.tsv)
if [[ "$case_name" != FeCl2 ]]; then arguments+=(--max-iterations "$cap"); fi
hostname > node.txt
module list > modules.txt 2>&1
env | sort > environment.txt
sha256sum "$binary" input.xmi *.gbs > hashes.txt
printf '%s\n' "${arguments[@]}" > arguments.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time "$binary" input.xmi "${arguments[@]}" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
# Native exit 2 is expected at the deliberate short-test iteration cap.
if [[ "$case_name" != FeCl2 && "$status" == 2 ]]; then exit 0; fi
exit "$status"
