#!/usr/bin/env bash
#SBATCH --job-name=ciah-metric-build
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=ciah-metric-build-%j.out
#SBATCH --error=ciah-metric-build-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
artifact="$SLURM_SUBMIT_DIR/ciah-metric-build-$SLURM_JOB_ID"
experiment="$SLURM_SUBMIT_DIR/ciah-metric.patch"
mkdir -p "$artifact"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
# The execution mirror is private to these build jobs. Restore its production
# source even on failure; this experiment never modifies the tracked repository.
patch --dry-run --batch -p1 -d "$source_dir" < "$experiment"
patch --batch -p1 -d "$source_dir" < "$experiment"
trap 'patch --batch -R -p1 -d "$source_dir" < "$experiment"' EXIT
cmake --build "$source_dir/build" --target xmvb -j8 > "$artifact/build.log" 2>&1
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/"
cp "$source_dir/build/src/CMakeFiles/xmvb.dir/flags.make" "$artifact/build-flags.txt"
cp "$experiment" "$artifact/"
cp "$SLURM_SUBMIT_DIR/ciah-positive-source.txt" "$artifact/source.txt"
sha256sum "$artifact/xmvb-cpp.exe" > "$artifact/binary.sha256"
