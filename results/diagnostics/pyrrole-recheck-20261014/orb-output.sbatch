#!/usr/bin/env bash
#SBATCH --job-name=orb-output-test
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=00:30:00
#SBATCH --output=orb-output-%j.out
#SBATCH --error=orb-output-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
run="$SLURM_SUBMIT_DIR/orb-output-$SLURM_JOB_ID"
mkdir -p "$run"
cd "$run"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cmake -S "$source_dir" -B "$source_dir/build" >configure.log 2>&1
cmake --build "$source_dir/build" --target xmvb test_output_orbitals -j16 >build.log 2>&1
srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores "$source_dir/build/src/test_output_orbitals" >unit.out 2>unit.err
fixture=/home/guqqgroup/xiatao/xmvb-runs/shared-response-20261012/ci-ray-regression-1252666/F2-OEO
cp "$fixture/input.xmi" F2.xmi
cp "$fixture/"*.gbs .
export OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores "$source_dir/build/src/xmvb-cpp.exe" F2.xmi --optimizer-backend neo --eigensolver davidson >F2.out 2>F2.err
test -s F2.orb
if grep -q '\$GUS\|\$END' F2.orb; then exit 1; fi
sha256sum "$source_dir/build/src/xmvb-cpp.exe" F2.xmi F2.orb >hashes.txt
