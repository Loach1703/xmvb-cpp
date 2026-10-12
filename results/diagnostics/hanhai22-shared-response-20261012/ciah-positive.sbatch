#!/usr/bin/env bash
#SBATCH --job-name=ciah-positive
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=ciah-positive-%A_%a.out
#SBATCH --error=ciah-positive-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -n ${BUILD_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
cases=(F2 241 240 MnF2 FeCl2)
case="${cases[$SLURM_ARRAY_TASK_ID]}"
artifact="$root/${RUN_LABEL:-ciah-positive}-$SLURM_ARRAY_JOB_ID/$case/one-electron"
frozen="$root/${BUILD_LABEL:-ciah-positive-build}-$BUILD_JOB_ID"
binary="$frozen/xmvb-cpp.exe"
mkdir -p "$artifact"
cp "$root/ciah-preconditioners-1247579/$case/input.xmi" "$artifact/"
cp "$root/regression-1246756/241/cc-pvdz.gbs" "$artifact/"
cp "$root/regression-1246756/FeCl2/def2-tzvp.gbs" "$artifact/"
cp "$frozen/build-flags.txt" "$frozen/source.txt" "$artifact/"
cd "$artifact"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > execution-env.txt
sha256sum "$binary" input.xmi *.gbs > hashes.txt
printf '%s\n' "$binary input.xmi --optimizer-backend neo --orbital-preconditioner one-electron --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --neo-trace neo.tsv" > command.txt
set +e
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
 "$binary" input.xmi --optimizer-backend neo --orbital-preconditioner one-electron \
 --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --neo-trace neo.tsv > output.out 2> run.err
status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
