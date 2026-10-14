#!/usr/bin/env bash
#SBATCH --job-name=porphine-screened-singlet
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=15-00:00:00
#SBATCH --output=porphine-screened-converge-%j.out
#SBATCH --error=porphine-screened-converge-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
run="$root/porphine-screened-converge-$SLURM_JOB_ID"
mkdir -p "$run"
cp "$root/inputs/singlet/input.xmi" "$root/inputs/singlet/cc-pvdz.gbs" "$run/"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
binary="$root/shell-screening-fix-1253607/xmvb-cpp.exe"
hostname > node.txt
sha256sum "$binary" input.xmi cc-pvdz.gbs > hashes.txt
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time "$binary" input.xmi \
  --optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal \
  --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 \
  --neo-trace neo.tsv > output.out 2> run.err
