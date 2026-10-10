#!/usr/bin/env bash
#SBATCH --job-name=mpi-service-check
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=15-00:00:00
#SBATCH --output=mpi-service-%j.out
#SBATCH --error=mpi-service-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06 openmpi/5.0.9
set -u
source_dir="$1"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=16 OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
for ranks in 1 2; do
    srun --exact --mpi=pmix_v5 --kill-on-bad-exit=1 --nodes="$ranks" \
        --ntasks="$ranks" --ntasks-per-node=1 --cpus-per-task=16 --cpu-bind=cores \
        "$source_dir/build/src/check_integral_service" "$source_dir/testdata/vbscf/F2.xmi"
done
# Missing CLI arguments exit without generating integrals; workers must shut down normally.
if srun --exact --mpi=pmix_v5 --nodes=2 --ntasks=2 --ntasks-per-node=1 \
    --cpus-per-task=16 --cpu-bind=cores "$source_dir/build/src/xmvb-cpp.exe"; then
    :
else
    status=$?
    [[ "$status" == 1 ]] || exit "$status"
fi
