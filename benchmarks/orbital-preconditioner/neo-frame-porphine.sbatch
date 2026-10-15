#!/usr/bin/env bash
#SBATCH --job-name=neo-frame-porphine
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=neo-frame-porphine-%A_%a.out
#SBATCH --error=neo-frame-porphine-%A_%a.err
set -eo pipefail
test -n "${SLURM_JOB_ID:-}"
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
frozen="$root/neo-frame-build-${1:?build job required}"
spins=(singlet triplet)
spin="${spins[$SLURM_ARRAY_TASK_ID]}"
run="$root/neo-frame-porphine-$SLURM_ARRAY_JOB_ID/$spin"
mkdir -p "$run"
cp "$root/inputs/$spin/input.xmi" "$root/inputs/$spin/cc-pvdz.gbs" "$run/"
cd "$run"
cp "$frozen/source.txt" "$frozen/flags.txt" .
hostname > node.txt
lscpu > cpu.txt
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export XMVB_NEO_TRACE=scalar XMVB_PHASE_TIMING=1
unset XMVB_NEO_CANDIDATE_AUDIT
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|XMVB_|SLURM_)' > execution-env.txt
binary="$frozen/xmvb-cpp.exe"
args=(input.xmi --optimizer-backend neo --eigensolver davidson
  --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3
  --energy-tolerance 1e-7 --exact-exchange-memory-mib 256
  --max-iterations 2 --neo-trace neo.tsv)
sha256sum "$binary" input.xmi cc-pvdz.gbs > hashes.txt
printf '%s ' "$binary" "${args[@]}" > command.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=verbose,cores \
  /usr/bin/time -v -o resource.time "$binary" "${args[@]}" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
# An unconverged two-step prefix normally exits 2. Inspect its termination
# reason; neither this script nor scheduler completion certifies convergence.
test "$status" = 0 || test "$status" = 2
