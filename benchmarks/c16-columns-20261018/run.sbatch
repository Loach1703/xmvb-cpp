#!/usr/bin/env bash
#SBATCH --job-name=C16-columns
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=00:20:00
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
export LC_ALL=C OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_MAX_ACTIVE_LEVELS=1
mkdir -p "$artifact/bin" "$artifact/native" "$artifact/phases"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
cmake --build "$source_dir/build-final" --target xmvb test_orthogonal_direct_ci \
  test_structure_tiles test_ci_transport test_selected_ci_transport test_coupled_structure -j16 \
  > "$artifact/build.log" 2>&1
ctest --test-dir "$source_dir/build-final" --output-on-failure \
  -R '^(orthogonal_direct_ci|structure_tiles|ci_transport|selected_ci_transport|coupled_structure)$' \
  > "$artifact/ctest.out" 2>&1
cmake --build "$source_dir/build-final" --target xmvb -j16 >> "$artifact/build.log" 2>&1
cp "$source_dir/build-final/src/xmvb-cpp.exe" "$artifact/bin/"
cp "$source_dir/build-final/CMakeCache.txt" "$artifact/"
cp "$SLURM_SUBMIT_DIR/source.sha256" "$artifact/"
cd "$source_dir"
sha256sum --check "$artifact/source.sha256" > "$artifact/source-check.txt"
for mode in native phases; do
  cp "$source_dir/benchmarks/c16-startup-20261018/C16_vb.xmi" "$artifact/$mode/"
  cp "$source_dir/basis/cc-pvdz.gbs" "$artifact/$mode/"
done
cd "$artifact"
sha256sum bin/xmvb-cpp.exe native/C16_vb.xmi native/cc-pvdz.gbs > hashes.txt
[[ ${C16_BUILD_ONLY:-0} != 1 ]] || exit 0
export OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export C16_ARTIFACT="$artifact"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
artifact="$C16_ARTIFACT"
grep Cpus_allowed_list /proc/self/status > "$artifact/affinity.txt"
for variable in OMP_NUM_THREADS OMP_PLACES OMP_PROC_BIND OMP_MAX_ACTIVE_LEVELS \
    OPENBLAS_NUM_THREADS MKL_NUM_THREADS BLIS_NUM_THREADS; do
  printf '%s=%s\n' "$variable" "${!variable-}"
done > "$artifact/threads.txt"
for mode in native phases; do
  cd "$artifact/$mode"
  unset XMVB_PHASE_TIMING
  export XMVB_CIAH_TRACE=scalar
  [[ $mode == native ]] || export XMVB_PHASE_TIMING=1
  set +e
  /usr/bin/time -v -o resource.time timeout -k 10s 120s \
    stdbuf -oL -eL "$artifact/bin/xmvb-cpp.exe" C16_vb.xmi --optimizer-backend ciah \
    > output.out 2> run.err
  printf '%s\n' "$?" > run.exit
  set -e
done
RUN
