#!/usr/bin/env bash
#SBATCH --job-name=ciah-preconditioners
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=ciah-preconditioners-%A_%a.out
#SBATCH --error=ciah-preconditioners-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
frozen="$root/ciah-1247545"
binary="$frozen/xmvb-cpp.exe"
cases=(F2 241 240 MnF2 FeCl2)
case="${cases[$SLURM_ARRAY_TASK_ID]}"
artifact="$root/ciah-preconditioners-$SLURM_ARRAY_JOB_ID/$case"
mkdir -p "$artifact"
if [[ "$case" == 240 ]]; then
 cp "$source_dir/testdata/vbscf/240_vb.xmi" "$artifact/input.xmi"
else
 cp "$root/regression-1246756/$case/input.xmi" "$artifact/input.xmi"
fi
cp "$root/regression-1246756/241/cc-pvdz.gbs" "$artifact/"
cp "$root/regression-1246756/FeCl2/def2-tzvp.gbs" "$artifact/"
cd "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > execution-env.txt
cp "$frozen/build-flags.txt" "$artifact/"
printf '%s\n' a8d60c8c5ac6c840a09f0fed26fd71a6f6c74238 > implementation-revision.txt
sha256sum "$binary" input.xmi *.gbs > hashes.txt
# Keep both variants on the same allocation and physical cores. Different
# systems run concurrently as independent array tasks; there are no repeats.
for preconditioner in hessian-diagonal one-electron; do
 mkdir -p "$preconditioner"
 cp input.xmi *.gbs "$preconditioner/"
 pushd "$preconditioner" >/dev/null
 printf '%s\n' "$binary input.xmi --optimizer-backend neo --orbital-preconditioner $preconditioner --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --neo-trace neo.tsv" > command.txt
 set +e
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time "$binary" input.xmi \
  --optimizer-backend neo --orbital-preconditioner "$preconditioner" \
  --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
  --neo-trace neo.tsv > output.out 2> run.err
 exit_code=$?
 set -e
 printf '%s\n' "$exit_code" > run.exit
 popd >/dev/null
done
