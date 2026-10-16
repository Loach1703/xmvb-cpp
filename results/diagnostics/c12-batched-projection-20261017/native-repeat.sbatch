#!/usr/bin/env bash
#SBATCH --job-name=c12-native-check
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --nodelist=anode40
#SBATCH --mem=64G
#SBATCH --output=c12-native-check-%j.out
#SBATCH --error=c12-native-check-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
set -u
root="$SLURM_SUBMIT_DIR"
frozen="$root/c12-regression-1264816"
artifact="$root/c12-native-check-$SLURM_JOB_ID"
mkdir -p "$artifact/native"
cp "$frozen/"{CMakeCache.txt,binaries.sha256,sources.sha256} "$artifact/"
cp "$frozen/native/"{input.xmi,cc-pvdz.gbs} "$artifact/native/"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE
export C12_NATIVE_ROOT="$artifact" C12_NATIVE_BINARY="$frozen/bin/xmvb-cpp.exe"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
grep Cpus_allowed_list /proc/self/status > "$C12_NATIVE_ROOT/affinity.txt"
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_CPUS_PER_TASK=)' > "$C12_NATIVE_ROOT/environment.txt"
sha256sum "$C12_NATIVE_BINARY" > "$C12_NATIVE_ROOT/executed-binary.sha256"
cd "$C12_NATIVE_ROOT/native"
sha256sum input.xmi cc-pvdz.gbs > input.sha256
/usr/bin/time -v -o resource.time "$C12_NATIVE_BINARY" input.xmi --ciah-trace ciah.tsv > output.out 2> run.err
grep -q '^Termination reason *: ciah_projected_gradient_tolerance' output.out
RUN
