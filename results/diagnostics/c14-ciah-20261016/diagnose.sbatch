#!/usr/bin/env bash
#SBATCH --job-name=C14-metric
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=128G
#SBATCH --time=00:30:00
#SBATCH --output=C14-metric-%j.out
#SBATCH --error=C14-metric-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
set -u
source_dir=/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015
artifact="$SLURM_SUBMIT_DIR/diagnose-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
cmake --build "$source_dir/build" --target xmvb -j8 > "$artifact/build.log" 2>&1
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/xmvb-cpp.exe"
cp "$SLURM_SUBMIT_DIR/C14.xmi" "$artifact/C14.xmi"
ln -s "$source_dir/basis" "$artifact/basis"
cd "$artifact"
sha256sum xmvb-cpp.exe C14.xmi > hashes.sha256
export XMVB_PHASE_TIMING=1 XMVB_CIAH_TRACE=scalar
set +e
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
 /usr/bin/time -v -o resource.time ./xmvb-cpp.exe C14.xmi \
 --ciah-trace ciah.tsv --max-iterations 1 > C14.out 2> run.err
status=$?
printf '%s\n' "$status" > run.exit
exit 0
