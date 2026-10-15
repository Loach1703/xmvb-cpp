#!/usr/bin/env bash
#SBATCH --job-name=porphine-quintet
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=15-00:00:00
#SBATCH --output=quintet-%A_%a.out
#SBATCH --error=quintet-%A_%a.err
set -eo pipefail
. /etc/profile
module purge
root="$SLURM_SUBMIT_DIR"
reference=/home/guqqgroup/xiatao/xmvb-runs/shared-response-20261012
modes=(neo pyscf-ah)
mode="${modes[$SLURM_ARRAY_TASK_ID]}"
run="$root/run-$SLURM_ARRAY_JOB_ID/$mode"
mkdir -p "$run"
cd "$run"
cp "$root/input.xmi" "$reference/inputs/triplet/cc-pvdz.gbs" .
export LC_ALL=C OMP_NUM_THREADS=32 OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
hostname > node.txt
lscpu > cpu.txt
sha256sum input.xmi cc-pvdz.gbs > hashes.txt
status=0
if [[ "$mode" == neo ]]; then
  module load gcc/14.1.0
  export OMP_PLACES=cores OMP_PROC_BIND=close
  export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
  export XMVB_NEO_TRACE=scalar
  unset XMVB_NEO_KEYFRAME_AUDIT XMVB_NEO_ORBITAL_STEP_AUDIT XMVB_NEO_ORBITAL_PATH XMVB_NEO_CANDIDATE_AUDIT
  cp /home/guqqgroup/xiatao/xmvb-runs/fecl2-fromcas-20261014/run-1255567/xmvb-cpp.exe .
  printf '%s  %s\n' cff4f5adb192a44b775086c95bd10fadf679231c37bd78a34ea63efd8d2e0100 xmvb-cpp.exe | sha256sum -c -
  sha256sum xmvb-cpp.exe >> hashes.txt
  srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
    ./xmvb-cpp.exe input.xmi --optimizer-backend neo --eigensolver davidson \
    --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
    --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err || status=$?
else
  module load anaconda/3/22.05
  export PYTHONPATH="$reference/porphine-reference-audit-1247977/python" PYTHONNOUSERSITE=1
  export OMP_PROC_BIND=false
  unset OMP_PLACES GOMP_CPU_AFFINITY
  cp "$root/porphine-cross-restart.py" .
  cmp input.xmi "$reference/inputs/quintet/input.xmi"
  sha256sum porphine-cross-restart.py >> hashes.txt
  srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
    python3 -u porphine-cross-restart.py "$reference" default-ah converge original quintet > output.out 2> run.err || status=$?
fi
printf '%s\n' "$status" > run.exit
exit "$status"
