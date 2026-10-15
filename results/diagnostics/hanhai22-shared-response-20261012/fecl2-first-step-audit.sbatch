#!/usr/bin/env bash
#SBATCH --job-name=fecl2-first-step-audit
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=fecl2-first-step-%j.out
#SBATCH --error=fecl2-first-step-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
artifact="$root/fecl2-first-step-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
cmake --build "$source_dir/build" --target check_exact_ctx_hvp -j8 > "$artifact/build.log" 2>&1
cp "$source_dir/build/src/check_exact_ctx_hvp" "$artifact/"
cp "$root/regression-1246756/FeCl2/input.xmi" "$artifact/"
cp "$root/regression-1246756/FeCl2/def2-tzvp.gbs" "$artifact/"
cd "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
sha256sum check_exact_ctx_hvp input.xmi > hashes.txt
export XMVB_NEO_TRACE=1
set +e
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  "$root/build-1246684/control.exe" input.xmi --max-iterations 1 \
  --dump-trace-dir trace --neo-trace neo.tsv \
  --dump-final-orbital-value-table-bin first-step.bin > first-step.out 2> first-step.err
status=$?
set -e
printf '%s\n' "$status" > first-step.exit
unset XMVB_NEO_TRACE
for sector in inactive-active inactive-external active-external; do
  srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    ./check_exact_ctx_hvp input.xmi --eigensolver davidson \
    --orbital-sector "$sector" --audit-oeo-chart 1 \
    --coupled-finite-difference 1 --step 0.0001 \
    > "$sector.out" 2> "$sector.err"
done
for sector in inactive-active inactive-external active-external; do
  srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    ./check_exact_ctx_hvp input.xmi --eigensolver davidson \
    --orbital-value-table-bin first-step.bin \
    --orbital-sector "$sector" --audit-oeo-chart 1 \
    --coupled-finite-difference 1 --step 0.0001 \
    > "first-$sector.out" 2> "first-$sector.err"
done
