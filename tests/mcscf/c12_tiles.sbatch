#!/usr/bin/env bash
#SBATCH --job-name=c12-tiles
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=64G
#SBATCH --output=c12-tiles-%j.out
#SBATCH --error=c12-tiles-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -d ${1:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
set -u
source_dir="$SLURM_SUBMIT_DIR"
baseline=$(realpath "$1")
artifact="$source_dir/c12-tiles-$SLURM_JOB_ID"
mkdir -p "$artifact/new"
export LC_ALL=C OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_CPUS_PER_TASK=)' > "$artifact/environment.txt"
cp "$baseline/provenance.txt" "$artifact/baseline-provenance.txt"
sha256sum "$baseline/"{test_structure_tiles,xmvb-cpp.exe} > "$artifact/baseline.sha256"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
 bash -c 'grep Cpus_allowed_list /proc/self/status' > "$artifact/affinity.txt"
find "$source_dir/src" "$source_dir/tests/vbscf/unit" -type f \
 \( -name '*.cpp' -o -name '*.hpp' \) -print0 | sort -z | xargs -0 sha256sum > "$artifact/sources.sha256"
cmake -S "$source_dir" -B "$source_dir/build" > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build" --target xmvb test_structure_tiles \
 test_representation_overlap test_selected_metric test_ci_transport \
 test_selected_ci_transport test_coupled_structure test_keyframe \
 test_orthogonal_direct_ci test_ciah -j16 > "$artifact/build.log" 2>&1
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
cp "$source_dir/build/src/"{test_structure_tiles,xmvb-cpp.exe} "$artifact/new/"
sha256sum "$artifact/new/"{test_structure_tiles,xmvb-cpp.exe} > "$artifact/new.sha256"
srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores env OMP_NUM_THREADS=2 \
 ctest --test-dir "$source_dir/build" --output-on-failure \
 -R '^(structure_tiles|representation_overlap|selected_metric|ci_transport|selected_ci_transport|coupled_structure|keyframe|orthogonal_direct_ci|ciah)$' \
 > "$artifact/ctest.out" 2>&1
cp "$source_dir/build/Testing/Temporary/LastTest.log" "$artifact/test-details.log"
for mode in selected12 structure-response12 product-jet12 gram12; do
 for variant in before after; do
  executable="$baseline/test_structure_tiles"
  [[ $variant == before ]] || executable="$artifact/new/test_structure_tiles"
  srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
   /usr/bin/time -v -o "$artifact/$variant-$mode.time" \
   "$executable" --"$mode" 132 > "$artifact/$variant-$mode.out" 2> "$artifact/$variant-$mode.err"
 done
done
# Independent numerical oracles are separate from paired native timing.
for mode in selected12-parity structure-response12-reference; do
 srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  "$artifact/new/test_structure_tiles" --"$mode" 132 \
  > "$artifact/$mode.out" 2> "$artifact/$mode.err"
done
# Real molecular end-to-end run: CAS(12,12), 132 covalent structures, exact ERI.
# Preserve the input guess and all tolerances; change only ISCF from 7 to 9.
for variant in before after; do
 case_dir="$artifact/$variant-loflea"
 mkdir -p "$case_dir"
 sed -E 's/^iscf=.*/iscf=9/I' "$source_dir/testdata/vbscf/csd_LOFLEA_vb.xmi" > "$case_dir/input.xmi"
 cp "$source_dir/basis/cc-pvdz.gbs" "$case_dir/"
 executable="$baseline/xmvb-cpp.exe"
 [[ $variant == before ]] || executable="$artifact/new/xmvb-cpp.exe"
 (
  cd "$case_dir"
  sha256sum input.xmi cc-pvdz.gbs > input.sha256
  srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
   /usr/bin/time -v -o resource.time "$executable" input.xmi --ciah-trace ciah.tsv \
   > output.out 2> run.err
 )
done
# Phase attribution is not used as a production wall-time measurement.
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
 env XMVB_PHASE_TIMING=1 "$artifact/new/test_structure_tiles" --gram12 132 \
 > "$artifact/gram-phases.out" 2> "$artifact/gram-phases.err"
