#!/usr/bin/env bash
#SBATCH --job-name=selected-ci
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=00:30:00
#SBATCH --output=selected-%j.out
#SBATCH --error=selected-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
deps=/home/guqqgroup/xiatao/software/xmvb-deps-20261012
artifact="$source_dir/selected-$SLURM_JOB_ID"
mkdir -p "$artifact"
export OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH="$deps/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
printenv OMP_NUM_THREADS OMP_PLACES OMP_PROC_BIND OMP_MAX_ACTIVE_LEVELS \
 OPENBLAS_NUM_THREADS MKL_NUM_THREADS > "$artifact/environment.txt"
sha256sum "$source_dir/src/mcscf/ci/product_paths.hpp" \
 "$source_dir/tests/vbscf/unit/test_structure_tiles.cpp" \
 "$source_dir/src/mcscf/ci/product_support.hpp" \
 "$source_dir/src/vbscf/structures/assembly/action.cpp" \
 "$source_dir/src/vbscf/structures/assembly/orthogonal_tiles.cpp" \
 "$source_dir/src/mcscf/ci/sigma_tiles.cpp" > "$artifact/sources.sha256"
if [[ ${1:-all} != benchmark ]]; then
cmake -S "$source_dir" -B "$source_dir/build" > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build" --target xmvb test_structure_tiles \
 test_orthogonal_direct_ci test_coupled_structure test_davidson test_neo \
 test_ciah test_ci_transport test_selected_ci_transport test_representation_overlap -j"$SLURM_CPUS_PER_TASK" \
 > "$artifact/build.log" 2>&1
sha256sum "$source_dir/build/src/xmvb-cpp.exe" \
 "$source_dir/build/src/test_structure_tiles" > "$artifact/binaries.sha256"
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
cd "$artifact"
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 ctest --test-dir "$source_dir/build" \
 -R '^(structure_tiles|orthogonal_direct_ci|coupled_structure|davidson|neo|ciah|ci_transport|selected_ci_transport|representation_overlap|neo_f2|eigensolver_equivalence_f2|eigensolver_equivalence_f2_oeo)$' \
 --verbose --output-on-failure > ctest.out 2>&1
if [[ ${1:-all} == validation ]]; then exit 0; fi
else
# Benchmark current source, including its narrow correctness suite, rather
# than trusting a binary left by an earlier job in the shared build directory.
cmake --build "$source_dir/build" --target test_structure_tiles -j8 > "$artifact/build.log" 2>&1
srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores \
 ctest --test-dir "$source_dir/build" -R '^structure_tiles$' \
 --output-on-failure > "$artifact/ctest.out" 2>&1
fi
cd "$artifact"
sha256sum "$source_dir/build/src/test_structure_tiles" > benchmark.sha256
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|SLURM_CPUS_PER_TASK=)' > benchmark_environment.txt
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 bash -c 'grep Cpus_allowed_list /proc/self/status' > affinity.txt
for count in 1 3 20 132; do
 if [[ ${1:-all} == paths ]]; then
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores /usr/bin/time -v -o "$count-paths.time" \
   "$source_dir/build/src/test_structure_tiles" --paths12 "$count" \
   > "$count-paths.out" 2> "$count-paths.err"
  continue
 fi
 srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores /usr/bin/time -v -o "$count.time" \
  "$source_dir/build/src/test_structure_tiles" --selected12 "$count" \
  > "$count.out" 2> "$count.err"
 srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores /usr/bin/time -v -o "$count-parity.time" \
  "$source_dir/build/src/test_structure_tiles" --selected12-parity "$count" \
  > "$count-parity.out" 2> "$count-parity.err"
done
