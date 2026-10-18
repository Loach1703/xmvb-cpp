#!/usr/bin/env bash
#SBATCH --job-name=C16-budgets
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=00:40:00
#SBATCH --output=budgets-%j.out
#SBATCH --error=budgets-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0
set -u
frozen=$(realpath -e "$1")
artifact="$SLURM_SUBMIT_DIR/budgets-$SLURM_JOB_ID"
bundle=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle
export LD_LIBRARY_PATH="$bundle/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_MAX_ACTIVE_LEVELS=1 OMP_PLACES=cores OMP_PROC_BIND=close
mkdir -p "$artifact/bin"
cp "$frozen/bin/xmvb-cpp.exe" "$artifact/bin/"
cp "$frozen/source.sha256" "$frozen/CMakeCache.txt" "$artifact/"
cp "$SLURM_SUBMIT_DIR/run_ci.py" "$artifact/"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
sha256sum "$artifact/bin/xmvb-cpp.exe" > "$artifact/binary.sha256"
export C16_ARTIFACT="$artifact"
export C16_INPUT="$frozen/native"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
artifact="$C16_ARTIFACT"
grep Cpus_allowed_list /proc/self/status > "$artifact/affinity.txt"
env | sort > "$artifact/environment.txt"
export XMVB_CIAH_TRACE=scalar
unset XMVB_PHASE_TIMING
# Same new executable, same node, sequential resource tradeoff measurements.
# Never rerun a historical executable. Each capture stops at initial CI solve.
for mib in 4096 1024 256; do
  mkdir -p "$artifact/mib-$mib"
  cp "$C16_INPUT/C16_vb.xmi" "$C16_INPUT/cc-pvdz.gbs" "$artifact/mib-$mib/"
  cd "$artifact/mib-$mib"
  python3 "$artifact/run_ci.py" --seconds 600 --workspace-mib "$mib" \
    "$artifact/bin/xmvb-cpp.exe" C16_vb.xmi > driver.out 2> driver.err
done
mkdir -p "$artifact/phases-4096"
cp "$C16_INPUT/C16_vb.xmi" "$C16_INPUT/cc-pvdz.gbs" "$artifact/phases-4096/"
cd "$artifact/phases-4096"
export XMVB_PHASE_TIMING=1
python3 "$artifact/run_ci.py" --seconds 120 --workspace-mib 4096 \
  "$artifact/bin/xmvb-cpp.exe" C16_vb.xmi > driver.out 2> driver.err
RUN
