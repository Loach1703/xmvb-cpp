#!/usr/bin/env bash
#SBATCH --job-name=coarse-porphine
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=01:00:00
#SBATCH --output=coarse-porphine-%A_%a.out
#SBATCH --error=coarse-porphine-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -n ${XMVB_COARSE_BUILD_JOB:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
spins=(singlet triplet)
spin="${spins[$SLURM_ARRAY_TASK_ID]}"
artifact="$root/coarse-porphine-$SLURM_ARRAY_JOB_ID/$spin"
mkdir -p "$artifact"
cd "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > node.txt
lscpu > cpu.txt
cp "$root/inputs/$spin/input.xmi" .
cp "$root/inputs/$spin/cc-pvdz.gbs" .
sed -i -E 's/itmax=[0-9]+/itmax=4/g' input.xmi
cp "$root/fecl2-ritz-cost-$XMVB_COARSE_BUILD_JOB/xmvb-cpp.exe" .
sha256sum xmvb-cpp.exe input.xmi cc-pvdz.gbs > hashes.txt
export XMVB_NEO_TRACE=scalar
unset XMVB_NEO_CANDIDATE_AUDIT
status=0
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
 ./xmvb-cpp.exe input.xmi --optimizer-backend neo --eigensolver davidson \
 --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
 --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
[[ "$status" == 0 || "$status" == 2 ]]
