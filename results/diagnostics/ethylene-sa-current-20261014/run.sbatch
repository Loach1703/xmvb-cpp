#!/usr/bin/env bash
#SBATCH --job-name=ethylene-sa-current
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=run-%j.out
#SBATCH --error=run-%j.err
set -eo pipefail
. /etc/profile
module purge
root="$SLURM_SUBMIT_DIR"
run="$root/run-$SLURM_JOB_ID"
mkdir -p "$run/cpp" "$run/pyscf"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
hostname > "$run/node.txt"
lscpu > "$run/cpu.txt"
cp "$root/input.xmi" "$root/aug-cc-pvdz.gbs" "$run/cpp/"
cp "$root/input.xmi" "$run/pyscf/"
cp "$root/reference.py" "$root/FeCl2_pyscf.py" "$run/pyscf/"
module load gcc/14.1.0
export OMP_PLACES=cores OMP_PROC_BIND=close
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
cd "$run/cpp"
cp /home/guqqgroup/xiatao/xmvb-runs/fecl2-fromcas-20261014/run-1255567/xmvb-cpp.exe .
printf '%s  %s\n' cff4f5adb192a44b775086c95bd10fadf679231c37bd78a34ea63efd8d2e0100 xmvb-cpp.exe | sha256sum -c -
sha256sum input.xmi xmvb-cpp.exe aug-cc-pvdz.gbs > hashes.txt
export XMVB_NEO_TRACE=scalar XMVB_PHASE_TIMING=1
unset XMVB_NEO_KEYFRAME_AUDIT XMVB_NEO_ORBITAL_STEP_AUDIT XMVB_NEO_ORBITAL_PATH XMVB_NEO_CANDIDATE_AUDIT
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
  ./xmvb-cpp.exe input.xmi --optimizer-backend neo --eigensolver davidson \
  --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
  --neo-trace neo.tsv > output.out 2> run.err
cd "$run/pyscf"
module purge
module load anaconda/3/22.05
export PYTHONPATH=/home/guqqgroup/xiatao/xmvb-runs/shared-response-20261012/porphine-reference-audit-1247977/python
export PYTHONNOUSERSITE=1 OMP_PROC_BIND=false
unset OMP_PLACES GOMP_CPU_AFFINITY XMVB_PHASE_TIMING XMVB_NEO_TRACE
sha256sum input.xmi reference.py FeCl2_pyscf.py > hashes.txt
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
  python3 -u reference.py input.xmi ah > output.out 2> run.err
