#!/usr/bin/env bash
#SBATCH --job-name=porphine-rotation
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=15-00:00:00
#SBATCH --output=porphine-rotation-%A_%a.out
#SBATCH --error=porphine-rotation-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -n ${XMVB_ROTATION_BUILD_JOB:-} ]] || exit 64
. /etc/profile
module purge
root="$SLURM_SUBMIT_DIR"
spins=(singlet singlet triplet triplet)
modes=(neo pyscf-ah neo pyscf-ah)
spin="${spins[$SLURM_ARRAY_TASK_ID]}"
mode="${modes[$SLURM_ARRAY_TASK_ID]}"
artifact="$root/porphine-rotation-$SLURM_ARRAY_JOB_ID/$spin/$mode"
mkdir -p "$artifact"
cd "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
hostname > node.txt
lscpu > cpu.txt
cp "$root/inputs/$spin/input.xmi" input.xmi
cp "$root/inputs/$spin/cc-pvdz.gbs" .
sha256sum input.xmi cc-pvdz.gbs > hashes.txt
if [[ "$mode" == neo ]]; then
 module load gcc/14.1.0
 export OMP_PLACES=cores OMP_PROC_BIND=close
 export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
 cp "$root/rotation-build-$XMVB_ROTATION_BUILD_JOB/xmvb-cpp.exe" .
 sha256sum xmvb-cpp.exe >> hashes.txt
 # No extra nonlinear keyframe audits; retain compact scalar diagnostics.
 export XMVB_NEO_TRACE=scalar
 unset XMVB_NEO_CANDIDATE_AUDIT
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
  ./xmvb-cpp.exe input.xmi --optimizer-backend neo --eigensolver davidson \
  --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
  --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err
else
 module load anaconda/3/22.05
 export PYTHONPATH="$root/porphine-reference-audit-1247977/python" PYTHONNOUSERSITE=1
 export OMP_PROC_BIND=false
 unset OMP_PLACES GOMP_CPU_AFFINITY
 cp "$root/porphine-cross-restart.py" .
 sha256sum porphine-cross-restart.py >> hashes.txt
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
  python3 -u porphine-cross-restart.py "$root" default-ah converge original "$spin" > output.out 2> run.err
fi
