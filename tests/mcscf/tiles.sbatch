#!/usr/bin/env bash
#SBATCH --job-name=orthogonal-tiles
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=00:30:00
#SBATCH --output=tiles-%j.out
#SBATCH --error=tiles-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
deps=/home/guqqgroup/xiatao/software/xmvb-deps-20261012
artifact="$source_dir/tiles-$SLURM_JOB_ID"
mkdir -p "$artifact"
export OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH="$deps/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
printenv OMP_NUM_THREADS OMP_PLACES OMP_PROC_BIND OMP_MAX_ACTIVE_LEVELS \
  OPENBLAS_NUM_THREADS MKL_NUM_THREADS > "$artifact/threads.txt"
cmake -S "$source_dir" -B "$source_dir/build" > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build" --target xmvb test_structure_tiles \
  test_orthogonal_direct_ci test_coupled_structure test_davidson test_neo -j16 \
  > "$artifact/build.log" 2>&1
sha256sum "$source_dir/build/src/xmvb-cpp.exe" \
  "$source_dir/build/src/test_structure_tiles" > "$artifact/binaries.sha256"
cd "$artifact"
export OMP_NUM_THREADS=2
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores ctest --test-dir "$source_dir/build" \
  -R '^(structure_tiles|orthogonal_direct_ci|coupled_structure|davidson|neo)$' \
  --output-on-failure > ctest.out 2>&1
export OMP_NUM_THREADS=32
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o kernel12.time \
  "$source_dir/build/src/test_structure_tiles" --benchmark12-no-reference \
  > kernel12.out 2> kernel12.err
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o kernel12-parity.time \
  "$source_dir/build/src/test_structure_tiles" --benchmark12 \
  > kernel12-parity.out 2> kernel12-parity.err
if [[ ${1:-molecule} == kernel ]]; then exit 0; fi
mkdir LOFLEA
cp "$source_dir/testdata/vbscf/csd_LOFLEA_vb.xmi" LOFLEA/input.xmi
cp "$source_dir/basis/cc-pvdz.gbs" LOFLEA/
cd LOFLEA
sha256sum input.xmi cc-pvdz.gbs > input.sha256
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -c \
  'grep Cpus_allowed_list /proc/self/status; scontrol show step "$SLURM_JOB_ID.$SLURM_STEP_ID"' \
  > affinity.txt
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
  "$source_dir/build/src/xmvb-cpp.exe" input.xmi \
  --optimizer-backend neo --eigensolver davidson \
  --orbital-preconditioner hessian-diagonal \
  --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
  --exact-exchange-memory-mib 256 --neo-trace neo.tsv > stdout.txt 2> stderr.txt
