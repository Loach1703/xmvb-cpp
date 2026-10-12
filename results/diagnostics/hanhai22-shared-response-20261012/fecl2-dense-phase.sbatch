#!/usr/bin/env bash
#SBATCH --job-name=fecl2-dense-phase
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=fecl2-dense-phase-%j.out
#SBATCH --error=fecl2-dense-phase-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -n ${XMVB_COARSE_BUILD_JOB:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
run="$root/fecl2-dense-phase-$SLURM_JOB_ID"
mkdir -p "$run"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
cp "$root/fecl2-ritz-cost-$XMVB_COARSE_BUILD_JOB/xmvb-cpp.exe" .
cp "$root/fecl2-metric-convergence-1248707/input.xmi" .
cp "$root/fecl2-metric-convergence-1248707/def2-tzvp.gbs" .
hostname > node.txt
lscpu > cpu.txt
sha256sum xmvb-cpp.exe input.xmi def2-tzvp.gbs > hashes.txt
export XMVB_NEO_TRACE=scalar
unset XMVB_NEO_CANDIDATE_AUDIT
status=0
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
 ./xmvb-cpp.exe input.xmi --optimizer-backend neo --eigensolver dense \
 --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
 --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
