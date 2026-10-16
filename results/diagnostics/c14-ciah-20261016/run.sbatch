#!/usr/bin/env bash
#SBATCH --job-name=C14-ciah
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=128G
#SBATCH --time=15-00:00:00
#SBATCH --output=C14-ciah-%j.out
#SBATCH --error=C14-ciah-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && ${SLURM_CPUS_PER_TASK:-0} == 32 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
set -u
source_dir=/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015
artifact="$SLURM_SUBMIT_DIR/run-$SLURM_JOB_ID"
mkdir -p "$artifact"
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/xmvb-cpp.exe"
cp "$SLURM_SUBMIT_DIR/C14.xmi" "$artifact/C14.xmi"
ln -s "$source_dir/basis" "$artifact/basis"
cd "$artifact"
printf '%s  xmvb-cpp.exe\n' \
  27966b5ea621c92392457a56ecd6e25d87167e3f31c56039333b64d963d70b57 \
  | sha256sum --check > binary-check.txt
sha256sum C14.xmi xmvb-cpp.exe basis/cc-pvdz.gbs > hashes.sha256
printf '%s\n' 48c145baea76357abb5105e38004f48dfdc8645b > production-revision.txt
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE
env | sort | sed -n '/^OMP_/p; /^OPENBLAS_/p; /^MKL_/p; /^BLIS_/p; /^SLURM_/p' > environment.txt
set +e
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -c '
  grep Cpus_allowed_list /proc/self/status > affinity.txt
  exec /usr/bin/time -v -o resource.time \
    ./xmvb-cpp.exe C14.xmi --ciah-trace ciah.tsv
' > C14.out 2> run.err
status=$?
printf '%s\n' "$status" > run.exit
set -e
[[ "$status" == 0 ]]
grep -q '^Termination reason *: ciah_projected_gradient_tolerance' C14.out
