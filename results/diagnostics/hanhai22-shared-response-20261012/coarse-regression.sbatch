#!/usr/bin/env bash
#SBATCH --job-name=coarse-regression
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=coarse-regression-%A_%a.out
#SBATCH --error=coarse-regression-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -n ${XMVB_COARSE_BUILD_JOB:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
cases=(F2 241 MnF2)
inputs=(ciah-preconditioners-1247579/F2 ciah-regression-1247549/241 ciah-preconditioners-1247579/MnF2)
case_name="${cases[$SLURM_ARRAY_TASK_ID]}"
input_dir="$root/${inputs[$SLURM_ARRAY_TASK_ID]}"
artifact="$root/coarse-regression-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$artifact"
cd "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > node.txt
lscpu > cpu.txt
cp "$input_dir/input.xmi" .
cp "$input_dir/"*.gbs .
sed -i -E 's/itmax=[0-9]+/itmax=1000/g' input.xmi
cp "$root/fecl2-ritz-cost-$XMVB_COARSE_BUILD_JOB/xmvb-cpp.exe" .
sha256sum xmvb-cpp.exe input.xmi ./*.gbs > hashes.txt
export XMVB_NEO_TRACE=scalar
unset XMVB_NEO_CANDIDATE_AUDIT
status=0
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
 ./xmvb-cpp.exe input.xmi --optimizer-backend neo --eigensolver davidson \
 --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
 --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
