#!/usr/bin/env bash
#SBATCH --job-name=c12-cost
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=64G
#SBATCH --output=c12-cost-%j.out
#SBATCH --error=c12-cost-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -d ${1:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
set -u
source_dir="$SLURM_SUBMIT_DIR"
bundle=$(realpath "$1")
artifact="$source_dir/c12-cost-$SLURM_JOB_ID"
mkdir -p "$artifact/bin"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export PATH="/home/guqqgroup/xiatao/software/gprofng-2.42/bin:$PATH"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/gprofng-2.42/lib:$LD_LIBRARY_PATH"
[[ ${2:-} == '' || ${2:-} == --capture-only ]] || exit 64
cp "$bundle/bin/"{test_structure_tiles,xmvb-cpp.exe} "$artifact/bin/"
cp "$bundle/"{CMakeCache.txt,sources.sha256} "$artifact/"
printf '%s\n' "$bundle" > "$artifact/bundle.txt"
export C12_COST_ROOT="$artifact"
export C12_COST_BUNDLE="$bundle" C12_CAPTURE_ONLY="${2:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
sha256sum "$artifact/bin/"* > "$artifact/binaries.sha256"
# All measurements use precisely the same CPU set, not merely the same node.
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
root="$C12_COST_ROOT"
bundle="$C12_COST_BUNDLE"
grep Cpus_allowed_list /proc/self/status > "$root/affinity.txt"
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_CPUS_PER_TASK=)' > "$root/environment.txt"
unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE
if [[ $C12_CAPTURE_ONLY != --capture-only ]]; then
for capture in native phases; do
 unset XMVB_PHASE_TIMING
 [[ $capture == native ]] || export XMVB_PHASE_TIMING=1
 for mode in selected12 structure-response12 product-jet12; do
  /usr/bin/time -v -o "$root/$capture-$mode.time" \
   "$root/bin/test_structure_tiles" --"$mode" 132 \
   > "$root/$capture-$mode.out" 2> "$root/$capture-$mode.err"
 done
done
mkdir -p "$root/loflea-phases"
cp "$bundle/native/"{input.xmi,cc-pvdz.gbs} "$root/loflea-phases/"
cd "$root/loflea-phases"
sha256sum input.xmi cc-pvdz.gbs > input.sha256
XMVB_PHASE_TIMING=1 XMVB_CIAH_TRACE=scalar /usr/bin/time -v -o resource.time \
 "$root/bin/xmvb-cpp.exe" input.xmi --ciah-trace ciah.tsv > output.out 2> run.err
fi
unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE
mkdir -p "$root/loflea-profile"
cp "$bundle/native/"{input.xmi,cc-pvdz.gbs} "$root/loflea-profile/"
cd "$root/loflea-profile"
sha256sum input.xmi cc-pvdz.gbs > input.sha256
profiler=''
perf=$(command -v perf || true)
if [[ -z $perf && -f /opt/intel/oneapi/vtune/latest/bin64/os-perf &&
      -x /opt/intel/oneapi/vtune/latest/bin64/os-perf ]]; then
 perf=/opt/intel/oneapi/vtune/latest/bin64/os-perf
fi
printf '%s\n' "$perf" > "$root/perf-path.txt"
if [[ -n $perf ]]; then
 if "$perf" stat -e cycles:u -- true > "$root/perf-probe.out" 2> "$root/perf-probe.err"; then
  profiler=perf; event=cycles:u
 elif "$perf" stat -e cpu-clock:u -- true > "$root/perf-software.out" 2> "$root/perf-software.err"; then
  profiler=perf; event=cpu-clock:u
 fi
fi
if [[ $profiler == perf ]]; then
 "$perf" record -o capture.perf -F 199 -e "$event" --call-graph fp -- \
  "$root/bin/xmvb-cpp.exe" input.xmi --ciah-trace ciah.tsv > output.out 2> run.err
 "$perf" report -i capture.perf --stdio --children > hotspots.txt 2> report.err
elif command -v gprofng > "$root/gprofng-path.txt" 2>&1; then
 profiler=gprofng
 gprofng --version > "$root/gprofng-version.txt" 2>&1
 /usr/bin/time -v -o collection.time gprofng collect app -p on -o capture.er \
  "$root/bin/xmvb-cpp.exe" input.xmi --ciah-trace ciah.tsv > output.out 2> run.err
 gprofng display text -header -statistics capture.er > quality.txt
 gprofng display text -limit 60 -functions -calltree capture.er > hotspots.txt
 gprofng display text -limit 60 -lines capture.er > lines.txt 2> lines.err
 gprofng display text -limit 64 -threads capture.er > threads.txt 2> threads.err
elif [[ -x /opt/intel/oneapi/vtune/latest/bin64/vtune ]]; then
 profiler=/opt/intel/oneapi/vtune/latest/bin64/vtune
 if "$profiler" -collect hotspots -knob sampling-mode=sw -result-dir capture.vtune -- \
     "$root/bin/xmvb-cpp.exe" input.xmi --ciah-trace ciah.tsv > output.out 2> run.err; then
  for report in summary hotspots callstacks; do
   "$profiler" -report "$report" -r capture.vtune > "$report.txt" 2> "$report.err"
  done
 else
  profiler=unavailable-vtune-rejected
 fi
else
 profiler=unavailable
fi
printf '%s\n' "$profiler" > "$root/profiler.txt"
RUN
