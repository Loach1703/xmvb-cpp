#!/usr/bin/env bash
#SBATCH --job-name=factor-gradient
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:20:00
#SBATCH --output=factor-gradient-%j.out
#SBATCH --error=factor-gradient-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -d ${1:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
source_bundle=$(realpath "$1")
artifact="$source_dir/factor-gradient-$SLURM_JOB_ID"
mkdir -p "$artifact/bundle"
export LC_ALL=C OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_CPUS_PER_TASK=)' > "$artifact/environment.txt"
cmake --build "$source_dir/build" --target xmvb test_orthogonal_direct_ci -j8 > "$artifact/build.log" 2>&1
cp "$source_dir/build/src/xmvb-cpp.exe" "$source_bundle/240_vb.xmi" \
 "$source_bundle/cc-pvdz.gbs" "$artifact/bundle/"
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
sha256sum "$artifact/bundle/xmvb-cpp.exe" > "$artifact/binary.sha256"
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 bash -c 'grep Cpus_allowed_list /proc/self/status' > "$artifact/affinity.txt"
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 "$source_dir/build/src/test_orthogonal_direct_ci" --benchmark-exchange \
 > "$artifact/sigma-exchange.csv" 2> "$artifact/sigma-exchange.err"
labels=(cov ion01 ion02 ion03 full)
keywords=('COV' 'ION(0:1)' 'ION(0:2)' 'ION(0:3)' 'FULL')
for index in 0 1 2 3 4; do
 mkdir -p "$artifact/${labels[index]}"
 cd "$artifact/${labels[index]}"
 sed -E -e "s/^str=.*/str=${keywords[index]}/I" -e 's/^iscf=.*/iscf=9/I' \
  "$artifact/bundle/240_vb.xmi" > input.xmi
 cp "$artifact/bundle/cc-pvdz.gbs" .
 sha256sum input.xmi cc-pvdz.gbs > input.sha256
 srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
  /usr/bin/time -v -o resource.time "$artifact/bundle/xmvb-cpp.exe" input.xmi \
  --optimizer-backend ciah --eigensolver davidson --orbital-preconditioner hessian-diagonal \
  --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --ciah-trace neo.tsv > output.out 2> run.err
 printf '0\n' > run.exit
done
