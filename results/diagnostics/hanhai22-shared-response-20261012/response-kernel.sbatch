#!/usr/bin/env bash
#SBATCH --job-name=response-kernel
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=response-kernel-%j.out
#SBATCH --error=response-kernel-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -n ${XMVB_RESPONSE_BUILD_JOB:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
artifact="$root/response-kernel-$SLURM_JOB_ID"
mkdir -p "$artifact"
cp /home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a/build/src/test_orthogonal_direct_ci "$artifact/kernel-test"
cp "$root/fecl2-ritz-cost-$XMVB_RESPONSE_BUILD_JOB/"{source.txt,source.patch,flags.txt} "$artifact/"
cd "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > execution-env.txt
sha256sum kernel-test > hashes.txt
printf '%s\n' 'srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores kernel-test --benchmark' > command.txt
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
  ./kernel-test --benchmark > kernels.csv 2> run.err
