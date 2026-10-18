#!/usr/bin/env bash
#SBATCH --job-name=CI-channels
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --array=0-4
#SBATCH --output=run-%A_%a.out
#SBATCH --error=run-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0
set -u
frozen=$(realpath -e "$1")
labels=(c12-cov c12-ion02 c12-full c14-cov c12-cov-64)
label=${labels[$SLURM_ARRAY_TASK_ID]}
input_case=$label
budget=4096
if [[ $label == c12-cov-64 ]]; then input_case=c12-cov; budget=64; fi
artifact="$SLURM_SUBMIT_DIR/run-$SLURM_ARRAY_JOB_ID/$label"
mkdir -p "$artifact/bin"
cp "$frozen/bin/xmvb-cpp.exe" "$artifact/bin/"
cp "$frozen/source.sha256" "$frozen/CMakeCache.txt" "$artifact/"
cp "$SLURM_SUBMIT_DIR/cases/$input_case/input.xmi" "$artifact/"
cp "$SLURM_SUBMIT_DIR/cc-pvdz.gbs" "$artifact/"
cd "$artifact"
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
export LD_LIBRARY_PATH=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle/lib:${LD_LIBRARY_PATH:-}
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE
if [[ ${CI_BENCH_TRACE:-} == scalar ]]; then export XMVB_CIAH_TRACE=scalar; fi
export CI_BENCH_ARTIFACT="$artifact" CI_BENCH_BUDGET="$budget"
env | sort > environment.txt
sha256sum bin/xmvb-cpp.exe input.xmi cc-pvdz.gbs > hashes.txt
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
cd "$CI_BENCH_ARTIFACT"
grep Cpus_allowed_list /proc/self/status > affinity.txt
set +e
/usr/bin/time -v -o resource.time ./bin/xmvb-cpp.exe input.xmi \
  --optimizer-backend ciah --ci-workspace-mib "$CI_BENCH_BUDGET" \
  --ciah-trace ciah.tsv > output.out 2> run.err
status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
RUN
