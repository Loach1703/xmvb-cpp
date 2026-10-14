#!/usr/bin/env bash
#SBATCH --job-name=fecl2-user-neo
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=run-%j.out
#SBATCH --error=run-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
mkdir -p "$root/run-$SLURM_JOB_ID"
cd "$root/run-$SLURM_JOB_ID"
cp "$root/input.xmi" "$root/def2-tzvp.gbs" .
binary=/home/guqqgroup/xiatao/xmvb-runs/fecl2-fromcas-20261014/run-1255567/xmvb-cpp.exe
printf '%s  %s\n' cff4f5adb192a44b775086c95bd10fadf679231c37bd78a34ea63efd8d2e0100 "$binary" | sha256sum -c -
cp "$binary" ./xmvb-cpp.exe
export OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname >node.txt
lscpu >cpu.txt
sha256sum xmvb-cpp.exe input.xmi def2-tzvp.gbs >hashes.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
 /usr/bin/time -v -o resource.time ./xmvb-cpp.exe input.xmi --neo-trace neo.tsv >output.out 2>run.err || status=$?
printf '%s\n' "$status" >run.exit
exit "$status"
