#!/usr/bin/env bash
#SBATCH --job-name=github-product12
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --sockets-per-node=1
#SBATCH --cores-per-socket=16
#SBATCH --threads-per-core=1
#SBATCH --mem=48G
#SBATCH --time=00:20:00
#SBATCH --output=github-product12-%j.out
#SBATCH --error=github-product12-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -d ${1:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
baseline=$(realpath "$1")
artifact="$source_dir/github-product12-$SLURM_JOB_ID"
build_root=${2:-$artifact}
dependency_root=/home/guqqgroup/xiatao/software/xmvb-deps-20261012
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=8 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="$dependency_root/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
env | sort > "$artifact/environment.txt"
cp "$baseline/source-commit.txt" "$artifact/baseline-commit.txt"
sha256sum "$source_dir/tests/mcscf/github_product12.cpp" > "$artifact/driver.sha256"
for variant in baseline current; do
 tree="$source_dir"
 baseline_mode=OFF
 if [[ $variant == baseline ]]; then tree="$baseline"; baseline_mode=ON; fi
 mkdir -p "$artifact/$variant"
 cmake -S "$tree" -B "$build_root/$variant/build" \
  -DBUILD_TESTING=OFF -DXMVB_BUILD_DEV_TARGETS=OFF -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG -g1 -fno-omit-frame-pointer' \
  -DXMVB_ENABLE_EIGEN_BLAS_BACKEND=ON \
  -DLAPACK_ROOT_DIR="$dependency_root" -DCINT_ROOT_DIR="$dependency_root" \
  -DEIGEN3_ROOT_DIR="$dependency_root" \
  -DCMAKE_PROJECT_INCLUDE="$source_dir/tests/mcscf/github_product12.cmake" \
  -DXMVB_BENCHMARK_BASELINE="$baseline_mode" > "$artifact/$variant/configure.log" 2>&1
 cmake --build "$build_root/$variant/build" --target benchmark_product12 -j8 \
  > "$artifact/$variant/build.log" 2>&1
 cp "$build_root/$variant/build/benchmark_product12" "$artifact/$variant/"
 cp "$build_root/$variant/build/CMakeCache.txt" "$artifact/$variant/"
 sha256sum "$artifact/$variant/benchmark_product12" > "$artifact/$variant/binary.sha256"
done
# Exactly one native run per variant/count. Both use the same eight LLC-local
# CPUs within this allocation; no historical executable is otherwise rerun.
for count in 1 3 20 132; do
 for variant in baseline current; do
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
   /usr/bin/time -v -o "$artifact/$variant/$count.time" \
   python3 "$source_dir/tests/mcscf/pin_ci.py" 8 "$artifact/cache-affinity.jsonl" \
   "$artifact/$variant/benchmark_product12" "$count" "$artifact/$variant/$count.bin" \
   > "$artifact/$variant/$count.out" 2> "$artifact/$variant/$count.err"
 done
done
# Attribute only the new implementation; never repeat the baseline for phases.
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 python3 "$source_dir/tests/mcscf/pin_ci.py" 8 "$artifact/cache-affinity.jsonl" \
 env XMVB_PHASE_TIMING=1 "$artifact/current/benchmark_product12" 132 \
 > "$artifact/current/phases.out" 2> "$artifact/current/phases.err"
