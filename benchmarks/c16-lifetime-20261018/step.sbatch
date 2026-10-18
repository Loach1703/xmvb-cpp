#!/usr/bin/env bash
#SBATCH --job-name=C16-one-step
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=00:35:00
#SBATCH --output=step-%j.out
#SBATCH --error=step-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0
set -u
frozen=$(realpath -e "$1")
artifact="$SLURM_SUBMIT_DIR/step-$SLURM_JOB_ID"
bundle=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle
export LD_LIBRARY_PATH="$bundle/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_MAX_ACTIVE_LEVELS=1 OMP_PLACES=cores OMP_PROC_BIND=close
mkdir -p "$artifact/bin"
cp "$frozen/bin/xmvb-cpp.exe" "$artifact/bin/"
cp "$frozen/source.sha256" "$frozen/CMakeCache.txt" "$artifact/"
cp "$frozen/native/C16_vb.xmi" "$frozen/native/cc-pvdz.gbs" "$artifact/"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
cd "$artifact"
sha256sum bin/xmvb-cpp.exe C16_vb.xmi cc-pvdz.gbs > hashes.txt
export XMVB_CIAH_TRACE=scalar XMVB_PHASE_TIMING=1
env | sort > environment.txt
export C16_ARTIFACT="$artifact"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
cd "$C16_ARTIFACT"
grep Cpus_allowed_list /proc/self/status > affinity.txt
set +e
timeout -k 10s 1800s /usr/bin/time -v -o resource.time \
  stdbuf -oL -eL ./bin/xmvb-cpp.exe C16_vb.xmi \
  --optimizer-backend ciah --ci-workspace-mib 4096 \
  --max-iterations 1 --ciah-trace ciah.tsv > output.out 2> run.err
status=$?
printf '%s\n' "$status" > exit-code.txt
# The one-step limit is intentional; completion is read from the accepted
# CIAH record, never inferred from the batch status or an iteration-limit exit.
RUN
