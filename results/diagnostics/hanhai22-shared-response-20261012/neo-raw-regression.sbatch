#!/usr/bin/env bash
#SBATCH --job-name=neo-raw-regression
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15-00:00:00
#SBATCH --output=neo-raw-regression-%A_%a.out
#SBATCH --error=neo-raw-regression-%A_%a.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
frozen="$root/neo-raw-build-1252034"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
cases=(FeCl2 MnF2 241 240 F2-OEO F2-SA2 F2-TN F2-block)
case_name="${cases[$SLURM_ARRAY_TASK_ID]}"
case "$case_name" in
  FeCl2|MnF2) reference="$root/inputs/$case_name" ;;
  241) reference="$root/joint-response-1251180/241" ;;
  240|F2-OEO) reference="$root/current-suite-1249844/$case_name" ;;
  *) reference="$root/inputs/F2" ;;
esac
run="$root/neo-raw-regression-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$run"
cp "$reference/input.xmi" "$run/"
shopt -s nullglob
assets=("$reference/"*.gbs)
if ((${#assets[@]})); then cp "${assets[@]}" "$run/"; fi
cd "$run"
backend=neo
case "$case_name" in
  F2-SA2) cp "$source_dir/testdata/vbscf/F2_SA2.xmi" input.xmi ;;
  F2-TN) backend=nonredundant_truncated_newton ;;
  F2-block) backend=block_lbfgs ;;
  F2-OEO) sed -i -E 's/ISCF=7/ISCF=9/g' input.xmi ;;
esac
hostname > node.txt
lscpu > cpu.txt
cp "$frozen/source.txt" "$frozen/flags.txt" .
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_NEO_TRACE XMVB_NEO_CANDIDATE_AUDIT XMVB_PHASE_TIMING
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|XMVB_|SLURM_)' > execution-env.txt
binary="$frozen/xmvb-cpp.exe"
sha256sum "$binary" input.xmi ./*.gbs > hashes.txt
args=(input.xmi --eigensolver davidson --gradient-tolerance 1e-3
  --energy-tolerance 1e-7 --exact-exchange-memory-mib 256)
if test "$case_name" != F2-OEO; then args+=(--optimizer-backend "$backend"); fi
if test "$backend" = neo; then
  args+=(--neo-trace neo.tsv)
  if test "$case_name" != F2-SA2; then args+=(--orbital-preconditioner hessian-diagonal); fi
fi
printf '%s ' "$binary" "${args[@]}" > command.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=verbose,cores \
  /usr/bin/time -v -o resource.time "$binary" "${args[@]}" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
