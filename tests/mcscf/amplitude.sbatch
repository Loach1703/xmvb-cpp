#!/usr/bin/env bash
#SBATCH --job-name=C12-amplitude
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=C12-amplitude-%j.out
#SBATCH --error=C12-amplitude-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
set -u
source_dir="$SLURM_SUBMIT_DIR"
qualified=$(realpath -e "$1")
artifact="$source_dir/c12-amplitude-$SLURM_JOB_ID"
mkdir -p "$artifact"
grep -q '100% tests passed, 0 tests failed out of 14' "$qualified/ctest.out"
sha256sum --check "$qualified/sources.sha256" > "$artifact/source-check.txt"
cp "$source_dir/build/src/test_structure_tiles" "$artifact/"
cp "$qualified/"{sources.sha256,CMakeCache.txt} "$artifact/"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export XMVB_PHASE_TIMING=1
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
sha256sum "$artifact/test_structure_tiles" > "$artifact/binary.sha256"
export C12_AMPLITUDE_ROOT="$artifact"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
cd "$C12_AMPLITUDE_ROOT"
grep Cpus_allowed_list /proc/self/status > affinity.txt
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > environment.txt
# Independent current-code full-carrier oracles, not old executable timings.
for mode in source-state-reference factor-direction-reference; do
 /usr/bin/time -v -o "$mode.time" ./test_structure_tiles --"$mode" 12 \
   > "$mode.out" 2> "$mode.err"
done
RUN
