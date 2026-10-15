#!/usr/bin/env bash
#SBATCH --job-name=raw-chart-suite
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15-00:00:00
#SBATCH --output=raw-chart-suite-%A_%a.out
#SBATCH --error=raw-chart-suite-%A_%a.err
set -eo pipefail
test -n "${SLURM_JOB_ID:-}"
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
frozen="$root/raw-chart-build-1251970"
cases=(F2-OEO 241 240 7975 QOCGUM CERRAS LOFLEA)
case_name="${cases[$SLURM_ARRAY_TASK_ID]}"
case "$case_name" in
  F2-OEO|240|7975) reference="$root/current-suite-1249844/$case_name" ;;
  LOFLEA) reference="$root/joint-isolated-1251207/$case_name" ;;
  *) reference="$root/joint-response-1251180/$case_name" ;;
esac
run="$root/raw-chart-suite-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$run"
cp "$reference/input.xmi" "$run/"
shopt -s nullglob
assets=("$reference/"*.gbs)
if ((${#assets[@]})); then cp "${assets[@]}" "$run/"; fi
cd "$run"
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
cp "$frozen/source.txt" "$frozen/flags.txt" .
printf '%s\n' "$reference" > reference-path.txt
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_NEO_TRACE XMVB_NEO_CANDIDATE_AUDIT XMVB_PHASE_TIMING
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|XMVB_|SLURM_)' > execution-env.txt
binary="$frozen/xmvb-cpp.exe"
sha256sum "$binary" input.xmi ./*.gbs > hashes.txt
args=(input.xmi --optimizer-backend neo --eigensolver davidson
  --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3
  --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 --neo-trace neo.tsv)
printf '%s ' "$binary" "${args[@]}" > command.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=verbose,cores \
  /usr/bin/time -v -o resource.time "$binary" "${args[@]}" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
