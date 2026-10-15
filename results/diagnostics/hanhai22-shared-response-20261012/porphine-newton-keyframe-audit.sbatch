#!/usr/bin/env bash
#SBATCH --job-name=porphine-newton-keyframe
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=01:00:00
#SBATCH --output=porphine-newton-keyframe-%j.out
#SBATCH --error=porphine-newton-keyframe-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0 anaconda/3/22.05
root="$SLURM_SUBMIT_DIR"
run="$root/porphine-newton-keyframe-$SLURM_JOB_ID"
mkdir -p "$run"
cp "$root/inputs/singlet/input.xmi" "$root/inputs/singlet/cc-pvdz.gbs" "$run/"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_PROC_BIND=false
unset OMP_PLACES GOMP_CPU_AFFINITY
export PYTHONPATH="$root/porphine-reference-audit-1247977/python" PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > node.txt
sha256sum "$root/porphine-newton-keyframe-audit.py" "$root/porphine-cross-restart.py" \
  input.xmi cc-pvdz.gbs > hashes.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time python3 -u "$root/porphine-newton-keyframe-audit.py" \
  "$root" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
