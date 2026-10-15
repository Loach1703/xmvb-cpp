#!/usr/bin/env bash
#SBATCH --job-name=porphine-horizontal
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --output=porphine-horizontal-%j.out
#SBATCH --error=porphine-horizontal-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
run="$root/porphine-horizontal-$SLURM_JOB_ID/singlet"
mkdir -p "$run"
cp "$root/inputs/singlet/input.xmi" "$root/inputs/singlet/"*.gbs "$run/"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export XMVB_NEO_TRACE=scalar XMVB_NEO_ORBITAL_STEP_AUDIT=1
unset XMVB_NEO_KEYFRAME_AUDIT XMVB_NEO_ORBITAL_PATH
binary="$root/horizontal-build-${BUILD_JOB:?}/xmvb-cpp.exe"
arguments=(--optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal
  --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 --neo-trace neo.tsv)
hostname > node.txt
env | sort > environment.txt
sha256sum "$binary" input.xmi *.gbs > hashes.txt
printf '%s\n' "${arguments[@]}" > arguments.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time "$binary" input.xmi "${arguments[@]}" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
