#!/usr/bin/env bash
#SBATCH --job-name=ciah-regression
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=ciah-regression-%A_%a.out
#SBATCH --error=ciah-regression-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -n ${CIAH_BUILD_JOB:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
cases=(F2 241 MnF2 FeCl2)
case="${cases[$SLURM_ARRAY_TASK_ID]}"
artifact="$root/ciah-regression-$SLURM_ARRAY_JOB_ID/$case"
mkdir -p "$artifact"
cp "$root/regression-1246756/$case/input.xmi" "$artifact/"
cp "$root/regression-1246756/$case/"*.gbs "$artifact/"
cd "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > node.txt
lscpu > cpu.txt
binary="$root/ciah-$CIAH_BUILD_JOB/xmvb-cpp.exe"
sha256sum "$binary" input.xmi *.gbs > hashes.txt
preconditioner=auto
[[ "$case" != FeCl2 ]] || preconditioner=one-electron
set +e
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
 /usr/bin/time -v -o resource.time "$binary" input.xmi \
 --optimizer-backend neo --orbital-preconditioner "$preconditioner" \
 --neo-trace neo.tsv > output.out 2> run.err
printf '%s\n' "$?" > run.exit
