#!/usr/bin/env bash
#SBATCH --job-name=C14-gram
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=128G
#SBATCH --time=00:30:00
#SBATCH --output=C14-gram-%j.out
#SBATCH --error=C14-gram-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
set -u
source_dir=/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015
artifact="$SLURM_SUBMIT_DIR/point-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
cmake -S "$source_dir" -B "$source_dir/build" > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build" --target test_selected_metric_point -j8 > "$artifact/build.log" 2>&1
cp "$source_dir/build/src/test_selected_metric_point" "$artifact/diagnostic.exe"
cp "$SLURM_SUBMIT_DIR/C14.xmi" "$artifact/C14.xmi"
ln -s "$source_dir/basis" "$artifact/basis"
cd "$artifact"
sha256sum diagnostic.exe C14.xmi > hashes.sha256
export XMVB_CIAH_TRACE=scalar
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
 /usr/bin/time -v -o resource.time ./diagnostic.exe C14.xmi point \
 > point.out 2> point.err
