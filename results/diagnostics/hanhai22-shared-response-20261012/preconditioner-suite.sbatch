#!/usr/bin/env bash
#SBATCH --job-name=ciah-precond-suite
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-14
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=preconditioner-suite-%A_%a.out
#SBATCH --error=preconditioner-suite-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
frozen="$root/ciah-positive-build-1247629"
binary="$frozen/xmvb-cpp.exe"
row=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$root/preconditioner-suite.tsv")
IFS=$'\t' read -r case input_source orbital_type <<< "$row"
[[ -n "$case" && -f "$root/$input_source" ]] || exit 64
artifact="$root/preconditioner-suite-$SLURM_ARRAY_JOB_ID/$case"
mkdir -p "$artifact"
cp "$root/$input_source" "$artifact/input.xmi"
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
printf '%s\n' "$orbital_type" > orbital-type.txt
sha256sum "$binary" input.xmi *.gbs > hashes.txt
# Balance run order across cases; paired variants always use the same node,
# executable, input and allocated cores. Do not impose a macro-iteration cap.
variants=(hessian-diagonal one-electron)
if (( SLURM_ARRAY_TASK_ID % 2 )); then variants=(one-electron hessian-diagonal); fi
for variant in "${variants[@]}"; do
 mkdir -p "$variant"
 cp input.xmi *.gbs "$variant/"
 pushd "$variant" >/dev/null
 printf '%s\n' "$binary input.xmi --optimizer-backend neo --eigensolver davidson --orbital-preconditioner $variant --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --neo-trace neo.tsv" > command.txt
 set +e
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
  "$binary" input.xmi --optimizer-backend neo --eigensolver davidson \
  --orbital-preconditioner "$variant" --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
  --neo-trace neo.tsv > output.out 2> run.err
 status=$?
 set -e
 printf '%s\n' "$status" > run.exit
 popd >/dev/null
done
