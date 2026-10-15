#!/usr/bin/env bash
#SBATCH --job-name=neo-current-suite
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-14
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=15-00:00:00
#SBATCH --output=current-suite-%A_%a.out
#SBATCH --error=current-suite-%A_%a.err

set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
frozen="$root/fecl2-ritz-cost-1249820"
binary="$frozen/xmvb-cpp.exe"
expected=9f51d9c85e5fabe3dcfd11108ff2075e86b60bb883b40b55d0922a8c70e2ab01
[[ $(sha256sum "$binary" | cut -d' ' -f1) == "$expected" ]] || exit 65
row=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$root/preconditioner-suite.tsv")
IFS=$'\t' read -r case_name input_source orbital_type <<< "$row"
[[ -n "$case_name" && -f "$root/$input_source" ]] || exit 64
artifact="$root/current-suite-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$artifact"
cp "$root/$input_source" "$artifact/input.xmi"
cp "$root/regression-1246756/241/"*.gbs "$artifact/"
cp "$root/regression-1246756/FeCl2/"*.gbs "$artifact/"
cp "$frozen/flags.txt" "$artifact/"
cd "$artifact"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export XMVB_NEO_TRACE=scalar
unset XMVB_NEO_CANDIDATE_AUDIT
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > execution-env.txt
printf '%s\n' "$orbital_type" > orbital-type.txt
printf '%s\n' '096d2419: restored balancing algorithm; opt-in origin audit disabled' > source.txt
sha256sum "$binary" input.xmi *.gbs > hashes.txt
printf '%s\n' "$binary input.xmi --optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 --neo-trace neo.tsv" > command.txt
status=0
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
 "$binary" input.xmi --optimizer-backend neo --eigensolver davidson \
 --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
 --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
