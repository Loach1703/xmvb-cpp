#!/usr/bin/env bash
#SBATCH --job-name=240-space-phases
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=00:20:00
#SBATCH --output=240-phase-%j.out
#SBATCH --error=240-phase-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -d ${1:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
export AUDIT_BUNDLE=$(realpath "$1")
export AUDIT_ROOT="$SLURM_SUBMIT_DIR/240-phase-$SLURM_JOB_ID"
mkdir -p "$AUDIT_ROOT"
if [[ ${2:-} == build ]]; then
 module load cmake/3.30.5
 cmake --build "$SLURM_SUBMIT_DIR/build" --target xmvb test_ci_transport test_selected_ci_transport \
  -j8 > "$AUDIT_ROOT/build.log" 2>&1
 mkdir "$AUDIT_ROOT/bundle"
 cp "$AUDIT_BUNDLE/240_vb.xmi" "$AUDIT_BUNDLE/cc-pvdz.gbs" "$AUDIT_ROOT/bundle/"
 cp "$SLURM_SUBMIT_DIR/build/src/xmvb-cpp.exe" "$AUDIT_ROOT/bundle/"
 export AUDIT_BUNDLE="$AUDIT_ROOT/bundle"
fi
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
hostname > "$AUDIT_ROOT/node.txt"
lscpu > "$AUDIT_ROOT/cpu.txt"
sha256sum "$AUDIT_BUNDLE/xmvb-cpp.exe" > "$AUDIT_ROOT/binary.sha256"
if [[ ${2:-} == build ]]; then
 srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores env OMP_NUM_THREADS=2 \
  ctest --test-dir "$SLURM_SUBMIT_DIR/build" -R '^(ci_transport|selected_ci_transport)$' \
  --output-on-failure > "$AUDIT_ROOT/ctest.out" 2>&1
fi
# One step holds exactly the same CPU set for all spaces and both measurement modes.
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores bash -se <<'RUN'
grep Cpus_allowed_list /proc/self/status > "$AUDIT_ROOT/affinity.txt"
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_|XMVB_)' > "$AUDIT_ROOT/environment.txt"
for mode in native phases; do
 for space in ion03 full; do
  mkdir -p "$AUDIT_ROOT/$mode/$space"
  cd "$AUDIT_ROOT/$mode/$space"
  keyword='ION(0:3)'
  if [[ $space == full ]]; then keyword=FULL; fi
  sed -E -e "s/^str=.*/str=$keyword/I" -e 's/^iscf=.*/iscf=9/I' "$AUDIT_BUNDLE/240_vb.xmi" > input.xmi
  cp "$AUDIT_BUNDLE/cc-pvdz.gbs" .
  unset XMVB_PHASE_TIMING
  if [[ $mode == phases ]]; then export XMVB_PHASE_TIMING=1; fi
  /usr/bin/time -v -o resource.time "$AUDIT_BUNDLE/xmvb-cpp.exe" input.xmi \
   --optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal \
   --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --neo-trace neo.tsv > output.out 2> run.err
 done
done
unset XMVB_PHASE_TIMING
cd "$AUDIT_ROOT"
if command -v perf > perf-path.txt 2>&1; then
 perf stat -e cpu-clock:u -- true > perf-probe.out 2> perf-probe.err && available=1 || available=0
else
 available=0
fi
if [[ $available == 1 ]]; then
 cd "$AUDIT_ROOT/phases/ion03"
 perf record -o "$AUDIT_ROOT/perf.data" -F 199 -e cpu-clock:u --call-graph fp -- \
  "$AUDIT_BUNDLE/xmvb-cpp.exe" input.xmi --optimizer-backend neo --eigensolver davidson \
  --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
  > "$AUDIT_ROOT/capture.out" 2> "$AUDIT_ROOT/capture.err"
 perf report -i "$AUDIT_ROOT/perf.data" --stdio --children > "$AUDIT_ROOT/perf-report.txt"
else
 command -v gprofng > "$AUDIT_ROOT/gprofng-path.txt" 2>&1 || true
fi
RUN
