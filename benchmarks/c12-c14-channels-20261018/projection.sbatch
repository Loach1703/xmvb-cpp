#!/usr/bin/env bash
#SBATCH --job-name=selected-projection-build
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=00:30:00
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
cmake --build "$source_dir/build-final" --target xmvb test_product_frontier \
  test_selected_ci_transport test_ci_transport test_coupled_structure \
  test_orthogonal_direct_ci test_representation_overlap test_selected_metric \
  test_ciah test_keyframe test_output_normalization -j16 > "$artifact/build.log" 2>&1
ctest --test-dir "$source_dir/build-final" --output-on-failure --timeout 120 \
  -R '^(product_frontier|selected_ci_transport|ci_transport|coupled_structure|orthogonal_direct_ci|representation_overlap|selected_metric|ciah|keyframe|output_normalization|text_report_f2|state_average_f2)$' \
  > "$artifact/ctest.out" 2>&1
cp "$source_dir/build-final/src/xmvb-cpp.exe" "$artifact/bin/"
cp "$SLURM_SUBMIT_DIR/source.sha256" "$source_dir/build-final/CMakeCache.txt" "$artifact/"
cd "$source_dir"
sha256sum --check "$artifact/source.sha256" > "$artifact/source-check.txt"
sha256sum "$artifact/bin/xmvb-cpp.exe" > "$artifact/binary.sha256"
