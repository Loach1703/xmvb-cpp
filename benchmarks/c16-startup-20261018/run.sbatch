#!/usr/bin/env bash
#SBATCH --job-name=C16-startup
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=00:15:00
#SBATCH --output=run-%j.out
#SBATCH --error=run-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0 cmake/3.31.9
set -u
source_dir=$(realpath -e "$1")
artifact="$SLURM_SUBMIT_DIR/run-$SLURM_JOB_ID"
bundle=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle
export LD_LIBRARY_PATH="$bundle/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
mkdir -p "$artifact/bin" "$artifact/native" "$artifact/profile"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
# The superseded P-LBFGS submission was cancelled during archive creation.
# Remove only that target's interrupted, zero-byte build products.
(cd "$source_dir/build-final/src" && \
  cmake -P CMakeFiles/xmvb_linalg.dir/cmake_clean.cmake) > "$artifact/archive-repair.log" 2>&1
cmake --build "$source_dir/build-final" --target xmvb -j16 > "$artifact/build.log" 2>&1
cp "$source_dir/build-final/src/xmvb-cpp.exe" "$artifact/bin/"
cp "$source_dir/build-final/CMakeCache.txt" "$artifact/"
cp "$SLURM_SUBMIT_DIR/source.sha256" "$artifact/"
cd "$source_dir"
sha256sum --check "$artifact/source.sha256" > "$artifact/source-check.txt"
for mode in native profile; do
  cp "$SLURM_SUBMIT_DIR/C16_vb.xmi" "$artifact/$mode/"
  cp "$source_dir/basis/cc-pvdz.gbs" "$artifact/$mode/"
done
cd "$artifact"
sha256sum bin/xmvb-cpp.exe native/C16_vb.xmi native/cc-pvdz.gbs > hashes.txt
export C16_ARTIFACT="$artifact"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
artifact="$C16_ARTIFACT"
grep Cpus_allowed_list /proc/self/status > "$artifact/affinity.txt"
env | sort > "$artifact/environment.txt"
cd "$artifact/native"
unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE
set +e
/usr/bin/time -v -o resource.time timeout -k 10s 120s \
  stdbuf -oL -eL "$artifact/bin/xmvb-cpp.exe" C16_vb.xmi --optimizer-backend ciah > output.out 2> run.err
printf '%s\n' "$?" > run.exit
set -e
cd "$artifact/profile"
export XMVB_PHASE_TIMING=1 XMVB_CIAH_TRACE=scalar
profiler=$(command -v gprofng || true)
if [[ -n $profiler ]]; then
  "$profiler" --version > version.txt 2>&1
  set +e
  timeout -k 10s 75s "$profiler" collect app -p on -o capture.er \
    stdbuf -oL -eL "$artifact/bin/xmvb-cpp.exe" C16_vb.xmi --optimizer-backend ciah > output.out 2> run.err
  printf '%s\n' "$?" > run.exit
  set -e
  "$profiler" display text -header -statistics capture.er > quality.txt 2> quality.err || true
  "$profiler" display text -limit 60 -functions -calltree capture.er > hotspots.txt 2> report.err || true
elif [[ -x /opt/intel/oneapi/vtune/latest/bin64/vtune ]]; then
  profiler=/opt/intel/oneapi/vtune/latest/bin64/vtune
  set +e
  timeout -k 15s 75s "$profiler" -collect hotspots -knob sampling-mode=sw \
    -result-dir capture.vtune -- "$artifact/bin/xmvb-cpp.exe" C16_vb.xmi --optimizer-backend ciah > output.out 2> run.err
  printf '%s\n' "$?" > run.exit
  set -e
  "$profiler" -finalize -r capture.vtune > finalize.out 2> finalize.err || true
  "$profiler" -report hotspots -r capture.vtune > hotspots.txt 2> report.err || true
  "$profiler" -report summary -r capture.vtune > quality.txt 2> quality.err || true
else
  profiler=unavailable
fi
printf '%s\n' "$profiler" > profiler.txt
RUN
