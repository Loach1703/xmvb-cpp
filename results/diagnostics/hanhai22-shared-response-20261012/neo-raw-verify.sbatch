#!/usr/bin/env bash
#SBATCH --job-name=neo-raw-verify
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=00:15:00
#SBATCH --output=neo-raw-verify-%j.out
#SBATCH --error=neo-raw-verify-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
artifact="$root/neo-raw-verify-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
cp "$root/neo-raw-quotient-test.cpp" "$source_dir/tests/vbscf/unit/test_sparse_orbital_quotient.cpp"
cp "$root/neo-raw-F2-PYSCF-LOW.xmi" "$source_dir/testdata/vbscf/F2_OEO_PYSCF_LOW.xmi"
sha256sum "$source_dir/tests/vbscf/unit/test_sparse_orbital_quotient.cpp" > "$artifact/source.sha256"
cmake --build "$source_dir/build" --target test_sparse_orbital_quotient -j16 > "$artifact/build.log" 2>&1
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  "$source_dir/build/src/test_sparse_orbital_quotient" > "$artifact/quotient.out" 2> "$artifact/quotient.err"
cd "$source_dir/build"
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  ctest -R '^neo_pyscf_reference_f2$' --output-on-failure > "$artifact/ctest.out" 2> "$artifact/ctest.err"
