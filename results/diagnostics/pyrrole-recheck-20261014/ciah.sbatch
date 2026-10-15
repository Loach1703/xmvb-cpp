#!/usr/bin/env bash
#SBATCH --job-name=sa-ciah
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --output=ciah-%A_%a.out
#SBATCH --error=ciah-%A_%a.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
cases=(F2-SA-HAO F2-SA-OEO oeo blo)
case_name="${cases[$SLURM_ARRAY_TASK_ID]}"
run="$root/ciah-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$run"
cp "$root/inputs/$case_name.xmi" "$run/input.xmi"
cp "$root/inputs/aug-cc-pvdz.gbs" "$root/inputs/cc-pvdz.gbs" "$run/"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export XMVB_NEO_TRACE=scalar
binary="$root/horizontal-build-${BUILD_JOB:?}/xmvb-cpp.exe"
sha256sum "$binary" input.xmi *.gbs > hashes.txt
hostname >node.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
 /usr/bin/time -v -o resource.time "$binary" input.xmi --optimizer-backend neo \
 --eigensolver davidson --orbital-preconditioner hessian-diagonal \
 --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 \
 --neo-trace neo.tsv >output.out 2>run.err || status=$?
printf '%s\n' "$status" > run.exit
if [[ "$case_name" == F2-SA-* && "$status" == 0 ]]; then
  srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    "$root/horizontal-build-$BUILD_JOB/check_exact_ctx_hvp" input.xmi \
    --eigensolver davidson --coupled-finite-difference 1 >derivatives.out 2>derivatives.err
fi
exit "$status"
