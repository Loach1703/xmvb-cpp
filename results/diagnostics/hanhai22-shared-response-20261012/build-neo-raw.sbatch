#!/usr/bin/env bash
#SBATCH --job-name=build-neo-raw
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=00:45:00
#SBATCH --output=build-neo-raw-%j.out
#SBATCH --error=build-neo-raw-%j.err
set -eo pipefail
test -n "${SLURM_JOB_ID:-}"
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
artifact="$root/neo-raw-build-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
sha256sum "$root/neo-raw-base.tar" "$root/neo-raw.patch" "$root/raw.cpp" "$root/raw.hpp" > "$artifact/source-inputs.sha256"
tar -xf "$root/neo-raw-base.tar" -C "$source_dir"
cd "$source_dir"
git apply --check "$root/neo-raw.patch"
git apply "$root/neo-raw.patch"
cp "$root/raw.cpp" "$root/raw.hpp" src/vbscf/orbitals/charts/
printf '%s\n' '2785651c plus production NEO unwhitened quotient, shared AO metric and reference-covector-only stopping. No diagnostic macros; TN/L-BFGS geometry unchanged. Unrelated RI edits excluded.' > "$artifact/source.txt"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
cmake --build build --target xmvb check_exact_ctx_hvp test_sparse_orbital_quotient test_ciah test_keyframe test_coupled_structure test_hessian_diagonal test_neo test_neo_response_solver test_neo_globalization test_block_inverse_bfgs test_oeo_normalization_pullback -j16 > "$artifact/build.log" 2>&1
for test_name in test_sparse_orbital_quotient test_ciah test_keyframe test_coupled_structure test_hessian_diagonal test_neo test_neo_response_solver test_neo_globalization test_block_inverse_bfgs test_oeo_normalization_pullback; do
  srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    "build/src/$test_name" > "$artifact/$test_name.out" 2> "$artifact/$test_name.err"
done
for case_name in FeCl2 MnF2; do
  cd "$root/inputs/$case_name"
  step=1e-4
  if test "$case_name" = MnF2; then step=1e-5; fi
  srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    "$source_dir/build/src/check_exact_ctx_hvp" input.xmi --eigensolver davidson \
    --step "$step" --coupled-finite-difference 1 > "$artifact/$case_name-hvp.out" 2> "$artifact/$case_name-hvp.err"
done
cd "$source_dir"
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  build/src/check_exact_ctx_hvp testdata/vbscf/F2_SA2.xmi --eigensolver davidson \
  --step 1e-4 --probe full --max-rel-error 1e-7 --coupled-finite-difference 1 > "$artifact/F2-SA2-hvp.out" 2> "$artifact/F2-SA2-hvp.err"
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/"
sha256sum "$artifact/xmvb-cpp.exe" > "$artifact/hash.txt"
grep -E 'CMAKE_BUILD_TYPE:|CMAKE_CXX_FLAGS(_RELEASE)?:' "$source_dir/build/CMakeCache.txt" > "$artifact/flags.txt"
