#!/usr/bin/env bash
#SBATCH --job-name=C14-current
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=64G
#SBATCH --output=C14-current-%j.out
#SBATCH --error=C14-current-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && ${SLURM_CPUS_PER_TASK:-0} == 32 && $# == 3 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
set -u
source_dir="$SLURM_SUBMIT_DIR"
input=$(realpath -e "$1")
qualified=$(realpath -e "$2")
revision="$3"
artifact="$source_dir/c14-current-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
printf '%s\n' "$revision" > "$artifact/revision.txt"
printf '%s\n' "$qualified" > "$artifact/qualification.txt"
# Reuse the already passed regression qualification, never an old SCF binary.
grep -q '100% tests passed' "$qualified/ctest.out"
sha256sum --check "$qualified/sources.sha256" > "$artifact/source-check.txt"
cp "$qualified/"{sources.sha256,ctest.out} "$artifact/"
cmake --build "$source_dir/build" --target xmvb -j16 > "$artifact/build.log" 2>&1
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/xmvb-cpp.exe"
cp "$input" "$artifact/C14.xmi"
cp "$source_dir/basis/cc-pvdz.gbs" "$artifact/cc-pvdz.gbs"
cd "$artifact"
sha256sum xmvb-cpp.exe C14.xmi cc-pvdz.gbs > hashes.sha256
env | sort | sed -n '/^OMP_/p; /^OPENBLAS_/p; /^MKL_/p; /^BLIS_/p; /^SLURM_/p' > environment.txt
set +e
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -c '
  grep Cpus_allowed_list /proc/self/status > affinity.txt
  exec /usr/bin/time -v -o resource.time ./xmvb-cpp.exe C14.xmi --ciah-trace ciah.tsv
' > C14.out 2> run.err
status=$?
printf '%s\n' "$status" > run.exit
set -e
[[ "$status" == 0 ]]
grep -q '^Termination reason *: ciah_projected_gradient_tolerance' C14.out
