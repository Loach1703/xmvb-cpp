#!/usr/bin/env bash
#SBATCH --job-name=horizontal-build
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=00:30:00
#SBATCH --output=horizontal-build-%j.out
#SBATCH --error=horizontal-build-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
run="$root/horizontal-build-$SLURM_JOB_ID"
mkdir -p "$run"
cd "$run"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
tests=(test_sparse_orbital_quotient test_selected_ci_transport test_ci_transport
  test_representation_overlap test_overlap_pullback test_oeo_normalization_pullback
  test_ciah test_coupled_structure test_normalized_orbital_curvature
  test_hessian_diagonal test_orbital_step_geometry test_shifted_metric_preconditioner)
cmake -S "$source_dir" -B "$source_dir/build" > configure.log 2>&1
cmake --build "$source_dir/build" --target xmvb "${tests[@]}" -j16 > build.log 2>&1
for test in "${tests[@]}"; do
  srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
    "$source_dir/build/src/$test" > "$test.out" 2> "$test.err"
done
cp "$source_dir/build/src/xmvb-cpp.exe" xmvb-cpp.exe
sha256sum xmvb-cpp.exe > hashes.txt
find "$source_dir/src/vbscf/orbitals/charts" -type f -exec sha256sum '{}' + > source-hashes.txt
sha256sum "$source_dir/src/vbscf/optimization/backends/ciah.cpp" \
  "$source_dir/tests/vbscf/unit/test_sparse_orbital_quotient.cpp" >> source-hashes.txt
