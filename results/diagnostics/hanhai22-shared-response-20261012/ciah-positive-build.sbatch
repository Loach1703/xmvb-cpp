#!/usr/bin/env bash
#SBATCH --job-name=ciah-positive-build
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=ciah-positive-build-%j.out
#SBATCH --error=ciah-positive-build-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
artifact="$SLURM_SUBMIT_DIR/ciah-positive-build-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
cmake -S "$source_dir" -B "$source_dir/build" > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build" --target xmvb test_ciah test_sparse_orbital_quotient test_shifted_metric_preconditioner -j8 > "$artifact/build.log" 2>&1
for test in test_ciah test_sparse_orbital_quotient test_shifted_metric_preconditioner; do
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores "$source_dir/build/src/$test" > "$artifact/$test.out" 2> "$artifact/$test.err"
done
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/"
cp "$source_dir/build/src/CMakeFiles/xmvb.dir/flags.make" "$artifact/build-flags.txt"
cp "$SLURM_SUBMIT_DIR/ciah-positive-source.txt" "$artifact/source.txt"
sha256sum "$artifact/xmvb-cpp.exe" > "$artifact/binary.sha256"
