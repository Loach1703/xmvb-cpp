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
if [[ ${2:-} == --build ]]; then
 # A cancelled compiler can leave a zero-byte object newer than its source.
 # Quarantine only those exact build artifacts before resuming make.
 while IFS= read -r -d '' object; do
  mv -- "$object" "$object.cancelled-$SLURM_JOB_ID"
 done < <(find "$source_dir/build/src/CMakeFiles" -type f -name '*.o' -empty -print0)
 cmake --build "$source_dir/build" --target xmvb test_structure_tiles \
  test_representation_overlap test_selected_ci_transport -j16 > "$artifact/build.log" 2>&1
 srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores env OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 \
  ctest --test-dir "$source_dir/build" --output-on-failure \
  -R '^(structure_tiles|representation_overlap|selected_ci_transport)$' > "$artifact/ctest.out" 2>&1
fi
cp "$bundle/"{test_structure_tiles,xmvb-cpp.exe} "$artifact/bin/"
export C12_COST_ROOT="$artifact"
export C12_COST_SOURCE="$source_dir"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
sha256sum "$artifact/bin/"* > "$artifact/binaries.sha256"
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
find "$source_dir/src" "$source_dir/tests/vbscf/unit" -type f \
 \( -name '*.cpp' -o -name '*.hpp' \) -print0 | sort -z | xargs -0 sha256sum > "$artifact/sources.sha256"
# All measurements use precisely the same CPU set, not merely the same node.
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
root="$C12_COST_ROOT"
source_dir="$C12_COST_SOURCE"
grep Cpus_allowed_list /proc/self/status > "$root/affinity.txt"
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_CPUS_PER_TASK=)' > "$root/environment.txt"
for capture in native phases; do
 unset XMVB_PHASE_TIMING
 [[ $capture == native ]] || export XMVB_PHASE_TIMING=1
 for mode in selected12 structure-response12 product-jet12; do
  /usr/bin/time -v -o "$root/$capture-$mode.time" \
   "$root/bin/test_structure_tiles" --"$mode" 132 \
   > "$root/$capture-$mode.out" 2> "$root/$capture-$mode.err"
 done
done
mkdir -p "$root/loflea-phases" "$root/loflea-profile"
for capture in phases profile; do
 cp "$source_dir/c12-tiles-1264541/after-loflea/"{input.xmi,cc-pvdz.gbs} "$root/loflea-$capture/"
done
cd "$root/loflea-phases"
sha256sum input.xmi cc-pvdz.gbs > input.sha256
XMVB_PHASE_TIMING=1 XMVB_CIAH_TRACE=scalar /usr/bin/time -v -o resource.time \
 "$root/bin/xmvb-cpp.exe" input.xmi --ciah-trace ciah.tsv > output.out 2> run.err
unset XMVB_PHASE_TIMING
cd "$root/loflea-profile"
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
elif command -v gprofng > "$root/gprofng-path.txt" 2>&1; then
 profiler=gprofng
 gprofng collect app -p on -o capture.er \
  "$root/bin/xmvb-cpp.exe" input.xmi --ciah-trace ciah.tsv > output.out 2> run.err
 gprofng display text -header -statistics capture.er > quality.txt
 gprofng display text -limit 60 -functions -calltree capture.er > hotspots.txt
else
 profiler=unavailable
fi
printf '%s\n' "$profiler" > "$root/profiler.txt"
RUN
