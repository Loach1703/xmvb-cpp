#!/usr/bin/env bash
#SBATCH --job-name=240-full32
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=00:15:00
#SBATCH --output=240-full32-%j.out
#SBATCH --error=240-full32-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -d ${1:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
export FULL32_INPUT=$(realpath "$1")
export FULL32_ROOT="$source_dir/full32-$SLURM_JOB_ID"
mkdir -p "$FULL32_ROOT"
export LC_ALL=C OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
# Reuse the current qualified executable, never a historical executable.
sha256sum -c "$source_dir/selected-transport-1262005/binary.sha256" > "$FULL32_ROOT/binary-check.txt"
sha256sum -c "$source_dir/selected-transport-1262005/sources.sha256" > "$FULL32_ROOT/source-check.txt"
cp "$source_dir/build/src/xmvb-cpp.exe" "$source_dir/build/CMakeCache.txt" "$FULL32_ROOT/"
cp "$source_dir/tests/mcscf/full32.provenance" "$FULL32_ROOT/"
sha256sum "$FULL32_ROOT/xmvb-cpp.exe" > "$FULL32_ROOT/binary.sha256"
sha256sum "$source_dir/tests/mcscf/full32.sbatch" > "$FULL32_ROOT/script.sha256"
ldd "$FULL32_ROOT/xmvb-cpp.exe" > "$FULL32_ROOT/libraries.txt"
hostname > "$FULL32_ROOT/node.txt"
lscpu > "$FULL32_ROOT/cpu.txt"
module list > "$FULL32_ROOT/modules.txt" 2>&1
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores bash -se <<'RUN'
grep Cpus_allowed_list /proc/self/status > "$FULL32_ROOT/affinity.txt"
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_|XMVB_)' > "$FULL32_ROOT/environment.txt"
for mode in native phases; do
 mkdir -p "$FULL32_ROOT/$mode/full"
 cd "$FULL32_ROOT/$mode/full"
 cp "$FULL32_INPUT/input.xmi" "$FULL32_INPUT/cc-pvdz.gbs" .
 sha256sum input.xmi cc-pvdz.gbs > input.sha256
 unset XMVB_PHASE_TIMING
 [[ $mode == phases ]] && export XMVB_PHASE_TIMING=1
 set +e
 /usr/bin/time -v -o resource.time "$FULL32_ROOT/xmvb-cpp.exe" input.xmi \
  --optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal \
  --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --neo-trace neo.tsv > output.out 2> run.err
 status=$?
 set -e
 printf '%s\n' "$status" > run.exit
 [[ $status == 0 ]] || exit "$status"
done
unset XMVB_PHASE_TIMING
mkdir -p "$FULL32_ROOT/profile"
cd "$FULL32_ROOT/profile"
cp "$FULL32_INPUT/input.xmi" "$FULL32_INPUT/cc-pvdz.gbs" .
available=0
if command -v perf > perf-path.txt 2>&1; then
 perf stat -e cpu-clock:u -- true > perf-probe.out 2> perf-probe.err && available=1 || true
fi
set +e
if [[ $available == 1 ]]; then
 perf record -o "$FULL32_ROOT/perf.data" -F 199 -e cpu-clock:u --call-graph fp -- \
  "$FULL32_ROOT/xmvb-cpp.exe" input.xmi --optimizer-backend neo --eigensolver davidson \
  --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
  > capture.out 2> capture.err
 status=$?
 printf '%s\n' "$status" > capture.exit
 perf report -i "$FULL32_ROOT/perf.data" --stdio --children > "$FULL32_ROOT/perf-report.txt" 2>&1
elif [[ -x /opt/intel/oneapi/vtune/latest/bin64/vtune ]]; then
 profiler=/opt/intel/oneapi/vtune/latest/bin64/vtune
 "$profiler" -version > profiler-version.txt 2>&1
 "$profiler" -collect hotspots -knob sampling-mode=sw -result-dir "$FULL32_ROOT/vtune" -- \
  "$FULL32_ROOT/xmvb-cpp.exe" input.xmi --optimizer-backend neo --eigensolver davidson \
  --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
  > capture.out 2> capture.err
 status=$?
 printf '%s\n' "$status" > capture.exit
 for report in summary hotspots callstacks; do
  "$profiler" -report "$report" -r "$FULL32_ROOT/vtune" \
   > "$FULL32_ROOT/vtune-$report.txt" 2> "$FULL32_ROOT/vtune-$report.err"
 done
else
 printf 'No permitted sampling profiler available.\n' > unavailable.txt
fi
exit 0
RUN
