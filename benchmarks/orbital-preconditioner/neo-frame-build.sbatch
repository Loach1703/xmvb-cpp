#!/usr/bin/env bash
#SBATCH --job-name=neo-frame-build
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=00:45:00
#SBATCH --output=neo-frame-build-%j.out
#SBATCH --error=neo-frame-build-%j.err
set -eo pipefail
test -n "${SLURM_JOB_ID:-}"
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
artifact="$root/neo-frame-build-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_NEO_TRACE XMVB_NEO_CANDIDATE_AUDIT XMVB_PHASE_TIMING
cd "$source_dir"
git apply --check "$root/neo-frame.patch"
git apply "$root/neo-frame.patch"
cp "$root/neo-frame-source.txt" "$root/neo-frame.patch" "$artifact/"
sha256sum "$root/neo-frame.patch" > "$artifact/source-inputs.sha256"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
cmake --build build --target xmvb test_sparse_orbital_quotient test_ciah test_keyframe \
  test_coupled_structure test_hessian_diagonal test_neo test_neo_response_solver test_neo_globalization \
  check_exact_ctx_hvp -j16 > "$artifact/build.log" 2>&1
for name in test_sparse_orbital_quotient test_ciah test_keyframe test_coupled_structure \
    test_hessian_diagonal test_neo test_neo_response_solver test_neo_globalization; do
  srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    "build/src/$name" > "$artifact/$name.out" 2> "$artifact/$name.err"
done
cd build
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  ctest -R '^neo_pyscf_reference_f2$' --output-on-failure > "$artifact/ctest.out" 2> "$artifact/ctest.err"
cp src/xmvb-cpp.exe "$artifact/"
cp "$artifact/neo-frame-source.txt" "$artifact/source.txt"
sha256sum "$artifact/xmvb-cpp.exe" > "$artifact/hash.txt"
grep -E 'CMAKE_BUILD_TYPE:|CMAKE_CXX_FLAGS(_RELEASE)?:' CMakeCache.txt > "$artifact/flags.txt"
