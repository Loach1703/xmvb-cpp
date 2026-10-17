#!/usr/bin/env bash
#SBATCH --job-name=cov12-step
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --sockets-per-node=1
#SBATCH --cores-per-socket=16
#SBATCH --threads-per-core=1
#SBATCH --mem=48G
#SBATCH --time=00:20:00
#SBATCH --output=cov12-step-%j.out
#SBATCH --error=cov12-step-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -n ${1:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
artifact="$source_dir/cov12-step-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=8 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
env | sort > "$artifact/environment.txt"
printf '%s\n' "$1" > "$artifact/revision.txt"
cmake --build "$source_dir/build" --target xmvb -j8 > "$artifact/build.log" 2>&1
cp "$source_dir/build/src/xmvb-cpp.exe" "$source_dir/build/CMakeCache.txt" "$artifact/"
sha256sum "$source_dir/src/mcscf/ci/"{sigma.hpp,sigma_tiles.cpp} \
 "$source_dir/src/vbscf/structures/assembly/product_jet.cpp" > "$artifact/sources.sha256"
# Preserve geometry, initial orbitals, COV, HAO, and packed-exact integrals.
sed -E 's/^iscf=.*/iscf=9/I' "$source_dir/testdata/vbscf/csd_LOFLEA_vb.xmi" > "$artifact/input.xmi"
cp "$source_dir/basis/cc-pvdz.gbs" "$artifact/"
cd "$artifact"
sha256sum input.xmi cc-pvdz.gbs xmvb-cpp.exe > hashes.sha256
set +e
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 /usr/bin/time -v -o resource.time \
 python3 "$source_dir/tests/mcscf/pin_ci.py" 8 "$artifact/affinity.jsonl" \
 "$artifact/xmvb-cpp.exe" input.xmi --optimizer-backend ciah --eigensolver davidson \
 --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 \
 --energy-tolerance 1e-7 --max-iterations 1 --ciah-trace neo.tsv > output.out 2> run.err
status=$?
set -e
printf '%s\n' "$status" > run.exit
# Reaching the requested one-step limit is not convergence; preserve the exit.
