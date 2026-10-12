#!/usr/bin/env bash
#SBATCH --job-name=porphine-low-branch
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=15-00:00:00
#SBATCH --output=porphine-low-branch-%A_%a.out
#SBATCH --error=porphine-low-branch-%A_%a.err
set -eo pipefail
. /etc/profile
module purge
root="$SLURM_SUBMIT_DIR"
modes=(neo pyscf-ah)
mode="${modes[$SLURM_ARRAY_TASK_ID]}"
artifact="$root/porphine-low-branch-$SLURM_ARRAY_JOB_ID/$mode"
mkdir -p "$artifact"
cd "$artifact"
export OMP_NUM_THREADS=32 OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LC_ALL=C
if [[ "$mode" == neo ]]; then
    module load gcc/14.1.0
    export OMP_PLACES=cores OMP_PROC_BIND=close
    export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
    # Freeze the already built and regressed final source binary.
    cp "$root/ciah-prefix-regression-1248081/xmvb-cpp.exe" .
    cp "$root/porphine-space-gradient-1247942/endpoint.xmi" restart.xmi
    cp "$root/inputs/triplet/cc-pvdz.gbs" .
    sha256sum xmvb-cpp.exe restart.xmi cc-pvdz.gbs > hashes.txt
    export XMVB_NEO_TRACE=scalar
    status=0
    srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
      ./xmvb-cpp.exe restart.xmi --optimizer-backend neo --eigensolver davidson \
      --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
      --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err || status=$?
    printf '%s\n' "$status" > run.exit
    exit "$status"
else
    module load anaconda/3/22.05
    export PYTHONPATH="$root/porphine-reference-audit-1247977/python" PYTHONNOUSERSITE=1
    export OMP_PROC_BIND=false
    unset OMP_PLACES GOMP_CPU_AFFINITY
    sha256sum "$root/porphine-cross-restart.py" > hashes.txt
    srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
      python3 -u "$root/porphine-cross-restart.py" "$root" default-ah converge > output.out 2> run.err
fi
