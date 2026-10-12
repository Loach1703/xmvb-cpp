#!/usr/bin/env bash
#SBATCH --job-name=rotation-build
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=rotation-build-%j.out
#SBATCH --error=rotation-build-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
artifact="$SLURM_SUBMIT_DIR/rotation-build-$SLURM_JOB_ID"
mkdir -p "$artifact"
export OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
cmake --build "$source_dir/build" --target xmvb test_ciah test_coupled_structure test_sparse_orbital_quotient -j8 > "$artifact/build.log" 2>&1
for target in test_ciah test_coupled_structure test_sparse_orbital_quotient; do
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores "$source_dir/build/src/$target" > "$artifact/$target.out" 2> "$artifact/$target.err"
done
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/"
sha256sum "$artifact/xmvb-cpp.exe" > "$artifact/hashes.txt"
