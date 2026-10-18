#!/usr/bin/env bash
#SBATCH --job-name=integral-transform
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=00:20:00
#SBATCH --output=qualify-%j.out
#SBATCH --error=qualify-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0 cmake/3.31.9
set -u
source_dir=$(realpath -e "$1")
artifact="$SLURM_SUBMIT_DIR/qualify-$SLURM_JOB_ID"
mkdir -p "$artifact/bin"
export LD_LIBRARY_PATH=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle/lib:${LD_LIBRARY_PATH:-}
export LC_ALL=C OMP_NUM_THREADS=2 OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
hostname > "$artifact/node.txt"
module list > "$artifact/modules.txt" 2>&1
cmake -S "$source_dir" -B "$source_dir/build-final" > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build-final" --target xmvb test_ri_ao_response \
  test_integral_transform test_ri_integral_transform test_orthogonal_ri_transform \
  test_ri_ao_h1e_hvp test_ri_active_two_electron_response test_ri_pair_contractions \
  check_exact_ctx_hvp test_product_frontier test_structure_tiles \
  test_selected_ci_transport test_ci_transport test_coupled_structure \
  test_orthogonal_direct_ci test_representation_overlap test_selected_metric \
  test_ciah test_keyframe test_hessian_diagonal test_output_normalization -j"$SLURM_CPUS_PER_TASK" > "$artifact/build.log" 2>&1
ctest --test-dir "$source_dir/build-final" --output-on-failure --timeout 120 \
  -R '^(integral_transform|ri_integral_transform|orthogonal_ri_transform|hessian_diagonal|ri_ao_response|ri_ao_h1e_hvp|ri_active_two_electron_response|ri_pair_contractions|ri_tnhvp_f2|product_frontier|selected_ci_transport|ci_transport|coupled_structure|orthogonal_direct_ci|structure_tiles|representation_overlap|selected_metric|ciah|keyframe|output_normalization|text_report_f2|state_average_f2)$' \
  > "$artifact/ctest.out" 2>&1
cp "$source_dir/build-final/src/xmvb-cpp.exe" "$artifact/bin/"
cp "$SLURM_SUBMIT_DIR/source.sha256" "$source_dir/build-final/CMakeCache.txt" "$artifact/"
cd "$source_dir"
sha256sum --check "$artifact/source.sha256" > "$artifact/source-check.txt"
sha256sum "$artifact/bin/xmvb-cpp.exe" > "$artifact/binary.sha256"
for mode in exact ri; do
for label in F2 FeCl2; do
  mkdir -p "$artifact/$mode-$label"
  cd "$artifact/$mode-$label"
  cp "$source_dir/basis/"{cc-pvdz.gbs,cc-pvdz-jkfit.gbs,def2-tzvp.gbs,def2-tzvp-jkfit.gbs} .
  "$source_dir/build-final/src/check_exact_ctx_hvp" "$source_dir/testdata/vbscf/$label.xmi" \
    --standard-two-electron-mode "$mode" --eigensolver dense --step 1e-4 \
    --max-rel-error 1e-5 --response-tolerance 1e-8 > hvp.out 2> hvp.err
done
done
