#!/usr/bin/env bash
#SBATCH --job-name=selected-metric
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:20:00
#SBATCH --output=selected-metric-%j.out
#SBATCH --error=selected-metric-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -d ${1:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
input_dir=$(realpath "$1")
artifact="$source_dir/selected-metric-$SLURM_JOB_ID"
mkdir -p "$artifact"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_CPUS_PER_TASK=)' > "$artifact/environment.txt"
cmake --build "$source_dir/build" --target xmvb -j8 > "$artifact/build.log" 2>&1
sha256sum "$source_dir/build/src/xmvb-cpp.exe" > "$artifact/binary.sha256"
sha256sum "$source_dir/src/vbscf/optimization/neo/ciah.cpp" \
 "$source_dir/src/vbscf/orbitals/charts/physical_metric.cpp" \
 "$source_dir/src/vbscf/optimization/trust_region/retraction.cpp" > "$artifact/sources.sha256"
cp "$input_dir/input.xmi" "$input_dir/cc-pvdz.gbs" "$artifact/"
cd "$artifact"
sha256sum input.xmi cc-pvdz.gbs > input.sha256
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 bash -c 'grep Cpus_allowed_list /proc/self/status' > affinity.txt
set +e
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 env XMVB_NEO_TRACE=scalar /usr/bin/time -v -o resource.time \
 "$source_dir/build/src/xmvb-cpp.exe" input.xmi \
 --optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal \
 --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --neo-trace neo.tsv \
 --dump-final-orbital-value-table-bin orbitals.bin > output.out 2> run.err
run_exit=$?
set -e
printf '%s\n' "$run_exit" > run.exit
