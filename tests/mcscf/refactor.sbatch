#!/usr/bin/env bash
#SBATCH --job-name=common-solvers
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:25:00
#SBATCH --output=common-solvers-%j.out
#SBATCH --error=common-solvers-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
artifact="$source_dir/common-solvers-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
find "$source_dir/src" "$source_dir/tests/mcscf/unit" \
 "$source_dir/tests/vbscf/unit" -type f \( -name '*.cpp' -o -name '*.hpp' \) \
 -print0 | sort -z | xargs -0 sha256sum > "$artifact/sources.sha256"
cmake --build "$source_dir/build" --target xmvb test_ci_factors test_casci \
 test_ciah test_neo test_neo_globalization test_minres test_block_inverse_bfgs \
 test_orthonormal_hvp_basis test_positive_ritz_secants test_shifted_metric_preconditioner \
 test_state_weights test_representation_overlap test_ci_transport test_selected_ci_transport \
 test_keyframe test_structure_tiles test_sparse_orbital_quotient test_coupled_structure \
 test_davidson test_eigen_response test_spectral_trust_region \
 audit_orbital_gradient benchmark_exact_ctx_hvp -j8 > "$artifact/build.log" 2>&1
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
sha256sum "$source_dir/build/src/xmvb-cpp.exe" "$source_dir/build/src/test_ciah" \
 "$source_dir/build/src/test_casci" > "$artifact/binaries.sha256"
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/xmvb-cpp.exe"
# These executables must contain no definitions from the VB model library.
for target in test_ciah test_neo test_minres test_casci; do
  nm -C "$source_dir/build/src/$target" > "$artifact/$target.symbols"
  if grep -E ' [TW] xmvb::vb::' "$artifact/$target.symbols"; then exit 65; fi
done
srun --exact --ntasks=1 --cpus-per-task=4 --cpu-bind=cores \
 ctest --test-dir "$source_dir/build" --output-on-failure \
 -R '^(ci_factors|casci|ciah|neo|neo_globalization|minres|block_inverse_bfgs|orthonormal_hvp_basis|positive_ritz_secants|shifted_metric_preconditioner|state_weights|representation_overlap|ci_transport|selected_ci_transport|keyframe|structure_tiles|sparse_orbital_quotient|coupled_structure|davidson|eigen_response|spectral_trust_region|neo_hessian_diagonal_f2_hao|neo_hessian_diagonal_f2_oeo)$' \
 > "$artifact/ctest.out" 2>&1
cp "$source_dir/build/Testing/Temporary/LastTest.log" "$artifact/test-details.log"
