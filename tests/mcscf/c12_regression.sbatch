#!/usr/bin/env bash
#SBATCH --job-name=c12-regression
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --nodelist=anode40
#SBATCH --mem=64G
#SBATCH --output=c12-regression-%j.out
#SBATCH --error=c12-regression-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
set -u
source_dir="$SLURM_SUBMIT_DIR"
artifact="$source_dir/c12-regression-$SLURM_JOB_ID"
mkdir -p "$artifact/bin"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
bundle=''
if [[ ${1:-} == --bundle ]]; then
 bundle=$(realpath "$2")
 # A queued run must not silently rebuild later source. Reuse the qualified
 # executable pair, original inputs and source manifest as one frozen bundle.
 cp "$bundle/"{CMakeCache.txt,sources.sha256} "$artifact/"
 cp "$bundle/bin/"{xmvb-cpp.exe,test_structure_tiles} "$artifact/bin/"
 printf '%s\n' "$bundle" > "$artifact/bundle.txt"
else
cmake -S "$source_dir" -B "$source_dir/build" > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build" --target xmvb test_structure_tiles \
 test_representation_overlap test_selected_metric test_ci_transport \
 test_selected_ci_transport test_coupled_structure test_keyframe \
 test_orthogonal_direct_ci test_ciah -j16 > "$artifact/build.log" 2>&1
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
find "$source_dir/src" "$source_dir/tests/vbscf/unit" "$source_dir/tests/linalg/unit" -type f \
 \( -name '*.cpp' -o -name '*.hpp' \) -print0 | sort -z | xargs -0 sha256sum > "$artifact/sources.sha256"
srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores env OMP_NUM_THREADS=2 \
 ctest --test-dir "$source_dir/build" --output-on-failure \
 -R '^(structure_tiles|representation_overlap|selected_metric|ci_transport|selected_ci_transport|coupled_structure|keyframe|orthogonal_direct_ci|ciah)$' \
 > "$artifact/ctest.out" 2>&1
cp "$source_dir/build/Testing/Temporary/LastTest.log" "$artifact/test-details.log"
if [[ ${1:-} == --tests-only ]]; then exit 0; fi
cp "$source_dir/build/src/"{xmvb-cpp.exe,test_structure_tiles} "$artifact/bin/"
fi
sha256sum "$artifact/bin/"* > "$artifact/binaries.sha256"
export C12_RUN_ROOT="$artifact" C12_RUN_SOURCE="$source_dir" C12_FROZEN_BUNDLE="$bundle"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
root="$C12_RUN_ROOT"
source_dir="$C12_RUN_SOURCE"
grep Cpus_allowed_list /proc/self/status > "$root/affinity.txt"
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_CPUS_PER_TASK=)' > "$root/environment.txt"
for mode in selected12-parity structure-response12-reference product-jet12; do
 "$root/bin/test_structure_tiles" --"$mode" 132 > "$root/$mode.out" 2> "$root/$mode.err"
done
"$root/bin/test_structure_tiles" --temporal12 > "$root/temporal12.out" 2> "$root/temporal12.err"
for mode in native loflea-phases; do
 mkdir -p "$root/$mode"
 cd "$root/$mode"
 if [[ -n $C12_FROZEN_BUNDLE ]]; then
  cp "$C12_FROZEN_BUNDLE/native/"{input.xmi,cc-pvdz.gbs} .
 else
  sed -E 's/^iscf=.*/iscf=9/I' "$source_dir/testdata/vbscf/csd_LOFLEA_vb.xmi" > input.xmi
  cp "$source_dir/basis/cc-pvdz.gbs" .
 fi
 sha256sum input.xmi cc-pvdz.gbs > input.sha256
 unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE
 if [[ $mode == loflea-phases ]]; then
  export XMVB_PHASE_TIMING=1 XMVB_CIAH_TRACE=scalar
 fi
 /usr/bin/time -v -o resource.time "$root/bin/xmvb-cpp.exe" input.xmi \
  --ciah-trace ciah.tsv > output.out 2> run.err
 grep -q '^Termination reason *: ciah_projected_gradient_tolerance' output.out
done
RUN
