#!/usr/bin/env bash
#SBATCH --job-name=porphine-master
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --output=porphine-master-%A_%a.out
#SBATCH --error=porphine-master-%A_%a.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
states=(singlet triplet)
state="${states[$SLURM_ARRAY_TASK_ID]}"
run="$root/porphine-master-$SLURM_ARRAY_JOB_ID/$state"
mkdir -p "$run"
cp "$root/inputs/$state/input.xmi" "$root/inputs/$state/"*.gbs "$run/"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export XMVB_NEO_TRACE=scalar
unset XMVB_NEO_KEYFRAME_AUDIT XMVB_NEO_ORBITAL_STEP_AUDIT XMVB_NEO_ORBITAL_PATH XMVB_NEO_CANDIDATE_AUDIT
binary="/home/guqqgroup/xiatao/xmvb-runs/pyrrole-recheck-20261014/horizontal-build-1255480/xmvb-cpp.exe"
expected="a5313100b12fd57e7d876e0d21c9333df5cc9d1a4d4dfca126c25b3f5da9b602"
printf '%s  %s\n' "$expected" "$binary" | sha256sum -c -
arguments=(--optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal
  --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 --neo-trace neo.tsv)
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
printf '%s\n' 'master 39201b81; build 1255480; local uncommitted RI edits excluded' > source.txt
printf '%s  %s\n' e7d58d653a5d99f6acee301e0096e582dedd7a9f47c22b4d2c2ace013ab351fb input.xmi | sha256sum -c -
env | sort > environment.txt
sha256sum "$binary" input.xmi *.gbs > hashes.txt
printf '%s\n' "${arguments[@]}" > arguments.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time "$binary" input.xmi "${arguments[@]}" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
