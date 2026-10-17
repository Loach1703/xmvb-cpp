#!/usr/bin/env bash
#SBATCH --job-name=joint-coupling
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --threads-per-core=1
#SBATCH --nodelist=anode40
#SBATCH --mem=16G
#SBATCH --output=joint-coupling-%j.out
#SBATCH --error=joint-coupling-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
export OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH=/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib
cmake --build build --target test_coupled_structure -j16
srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores \
 ctest --test-dir build --output-on-failure -R '^coupled_structure$'
sha256sum tests/vbscf/unit/test_coupled_structure.cpp build/src/test_coupled_structure
