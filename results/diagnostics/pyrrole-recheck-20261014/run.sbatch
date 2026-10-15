#!/usr/bin/env bash
#SBATCH --job-name=pyrrole-recheck
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --output=pyrrole-%A_%a.out
#SBATCH --error=pyrrole-%A_%a.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
types=(oeo blo oeo blo)
methods=(block-lbfgs block-lbfgs neo neo)
type="${types[$SLURM_ARRAY_TASK_ID]}"
method="${methods[$SLURM_ARRAY_TASK_ID]}"
run="$root/run-$SLURM_ARRAY_JOB_ID/$type-$method"
mkdir -p "$run"
cp "$root/inputs/$type.xmi" "$run/input.xmi"
cp "$root/inputs/aug-cc-pvdz.gbs" "$run/"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export XMVB_NEO_TRACE=scalar
binary=/home/guqqgroup/xiatao/xmvb-runs/shared-response-20261012/horizontal-build-1255164/xmvb-cpp.exe
printf '%s  %s\n' bc2003c898dddf977d76ac6bb6e4a5fd03f0700dd9957595c2eb73e830cfcfbb "$binary" | sha256sum -c -
arguments=(--eigensolver davidson --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256)
if [[ "$method" == neo ]]; then
  arguments+=(--optimizer-backend neo --orbital-preconditioner hessian-diagonal --neo-trace neo.tsv)
fi
hostname > node.txt
env | sort > environment.txt
sha256sum "$binary" input.xmi aug-cc-pvdz.gbs > hashes.txt
printf '%s\n' "${arguments[@]}" > arguments.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time "$binary" input.xmi "${arguments[@]}" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
