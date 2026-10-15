#!/usr/bin/env bash
#SBATCH --job-name=porphine-main-neo
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --exclusive
#SBATCH --time=15-00:00:00
#SBATCH --output=porphine-main-%A_%a.out
#SBATCH --error=porphine-main-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -n ${XMVB_PORPHINE_BUILD_JOB:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
spins=(singlet triplet)
spin="${spins[$SLURM_ARRAY_TASK_ID]}"
frozen="$root/porphine-main-build-$XMVB_PORPHINE_BUILD_JOB"
artifact="$root/porphine-main-$SLURM_ARRAY_JOB_ID/$spin"
mkdir -p "$artifact"
cp "$root/inputs/$spin/input.xmi" "$artifact/"
cp "$root/inputs/$spin/cc-pvdz.gbs" "$artifact/"
cp "$frozen/xmvb-cpp.exe" "$frozen/source.txt" "$frozen/flags.txt" "$artifact/"
cd "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_NEO_TRACE XMVB_NEO_CANDIDATE_AUDIT XMVB_PHASE_TIMING
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > execution-env.txt
sha256sum xmvb-cpp.exe input.xmi cc-pvdz.gbs > hashes.txt
printf '%s\n' './xmvb-cpp.exe input.xmi --optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 --neo-trace neo.tsv' > command.txt
status=0
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=verbose,cores /usr/bin/time -v -o resource.time \
 ./xmvb-cpp.exe input.xmi --optimizer-backend neo --eigensolver davidson \
 --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
 --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
