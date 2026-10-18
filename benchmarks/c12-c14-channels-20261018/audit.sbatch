#!/usr/bin/env bash
#SBATCH --job-name=CI-end-to-end-audit
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=00:30:00
#SBATCH --output=audit-%j.out
#SBATCH --error=audit-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 3 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0
set -u
frozen=$(realpath -e "$1")
label=$2
mode=$3
[[ $label == c12-cov || $label == c12-ion02 || $label == c12-full ]] || exit 64
[[ $mode == phase || $mode == profile ]] || exit 64
artifact="$SLURM_SUBMIT_DIR/audit-$SLURM_JOB_ID"
mkdir -p "$artifact/bin"
cp "$frozen/bin/xmvb-cpp.exe" "$artifact/bin/"
cp "$frozen/source.sha256" "$frozen/CMakeCache.txt" "$artifact/"
cp "$SLURM_SUBMIT_DIR/cases/$label/input.xmi" "$artifact/"
cp "$SLURM_SUBMIT_DIR/cc-pvdz.gbs" "$artifact/"
cd "$artifact"
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
export LD_LIBRARY_PATH=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle/lib:${LD_LIBRARY_PATH:-}
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export XMVB_PHASE_TIMING=1 XMVB_CIAH_TRACE=scalar
export CI_AUDIT_ARTIFACT="$artifact" CI_AUDIT_MODE="$mode"
env | sort > environment.txt
sha256sum bin/xmvb-cpp.exe input.xmi cc-pvdz.gbs > hashes.txt
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
cd "$CI_AUDIT_ARTIFACT"
grep Cpus_allowed_list /proc/self/status > affinity.txt
command=(./bin/xmvb-cpp.exe input.xmi --optimizer-backend ciah
         --ci-workspace-mib 4096 --ciah-trace ciah.tsv)
if [[ $CI_AUDIT_MODE == profile ]]; then
  gprofng --version > profiler-version.txt
  command=(gprofng collect app -p on -o capture.er "${command[@]}")
fi
set +e
/usr/bin/time -v -o resource.time "${command[@]}" > output.out 2> run.err
status=$?
printf '%s\n' "$status" > run.exit
if [[ $CI_AUDIT_MODE == profile ]]; then
  gprofng display text -header -statistics capture.er > quality.txt 2>&1
  gprofng display text -limit 60 -functions capture.er > functions.txt 2>&1
  gprofng display text -limit 100 -calltree capture.er > calltree.txt 2>&1
fi
exit "$status"
RUN
