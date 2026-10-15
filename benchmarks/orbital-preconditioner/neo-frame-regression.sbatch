#!/usr/bin/env bash
#SBATCH --job-name=neo-frame-regression
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=neo-frame-regression-%A_%a.out
#SBATCH --error=neo-frame-regression-%A_%a.err
set -eo pipefail
test -n "${SLURM_JOB_ID:-}"
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
frozen="$root/neo-frame-build-${1:?build job required}"
cases=(FeCl2 MnF2 241 240 F2-OEO F2-SA2 F2-TN F2-block)
case_name="${cases[$SLURM_ARRAY_TASK_ID]}"
case "$case_name" in
  FeCl2|MnF2) reference="$root/neo-raw-regression-1252070/$case_name" ;;
  241) reference="$root/neo-raw-regression-1252070/241" ;;
  240|F2-OEO|F2-SA2|F2-TN) reference="$root/neo-raw-regression-1252077/$case_name" ;;
  F2-block) reference="$root/neo-raw-regression-1252083/F2-block" ;;
esac
run="$root/neo-frame-regression-$SLURM_ARRAY_JOB_ID/$case_name"
mkdir -p "$run"
cp "$reference/input.xmi" "$run/"
shopt -s nullglob
assets=("$reference/"*.gbs)
if ((${#assets[@]})); then cp "${assets[@]}" "$run/"; fi
cd "$run"
hostname > node.txt
lscpu > cpu.txt
cp "$frozen/source.txt" "$frozen/flags.txt" .
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_NEO_TRACE XMVB_NEO_CANDIDATE_AUDIT XMVB_PHASE_TIMING
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|XMVB_|SLURM_)' > execution-env.txt
binary="$frozen/xmvb-cpp.exe"
args=(input.xmi --eigensolver davidson --gradient-tolerance 1e-3
  --energy-tolerance 1e-7 --exact-exchange-memory-mib 256)
case "$case_name" in
  F2-TN) args+=(--optimizer-backend nonredundant_truncated_newton) ;;
  F2-block) args+=(--optimizer-backend block_lbfgs) ;;
  *) args+=(--optimizer-backend neo --neo-trace neo.tsv)
     if test "$case_name" != F2-SA2; then args+=(--orbital-preconditioner hessian-diagonal); fi ;;
esac
sha256sum "$binary" input.xmi ./*.gbs > hashes.txt
printf '%s ' "$binary" "${args[@]}" > command.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=verbose,cores \
  /usr/bin/time -v -o resource.time "$binary" "${args[@]}" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
# A separate scalar capture counts remaining factorizations without new physics.
if test "$case_name" = FeCl2; then
  mkdir phases
  cp input.xmi ./*.gbs phases/
  cd phases
  export XMVB_PHASE_TIMING=1 XMVB_NEO_TRACE=scalar
  phase_status=0
  srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    "$binary" "${args[@]}" --max-iterations 2 > output.out 2> run.err || phase_status=$?
  printf '%s\n' "$phase_status" > run.exit
  test "$phase_status" = 0 || test "$phase_status" = 2
fi
exit "$status"
