#!/usr/bin/env bash
#SBATCH --job-name=porphine-path
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=01:00:00
#SBATCH --output=porphine-path-%A_%a.out
#SBATCH --error=porphine-path-%A_%a.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0 anaconda/3/22.05
root="$SLURM_SUBMIT_DIR"
tags=(neo-1 neo-4 pyscf-4 pyscf-newton-1)
tag="${tags[$SLURM_ARRAY_TASK_ID]}"
run="$root/porphine-path-$SLURM_ARRAY_JOB_ID/$tag"
mkdir -p "$run"
cp "$root/inputs/singlet/input.xmi" "$root/inputs/singlet/cc-pvdz.gbs" "$run/"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > node.txt
status=0
if [[ "$tag" == pyscf-* ]]; then
  export PYTHONPATH="$root/porphine-reference-audit-1247977/python" PYTHONNOUSERSITE=1
  export OMP_PROC_BIND=false
  unset OMP_PLACES GOMP_CPU_AFFINITY
  if [[ "$tag" == pyscf-4 ]]; then
    script="$root/porphine-path-pyscf.py"
  else
    script="$root/porphine-path-newton.py"
  fi
  sha256sum "$script" "$root/porphine-cross-restart.py" \
    input.xmi cc-pvdz.gbs > hashes.txt
  srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    /usr/bin/time -v -o resource.time python3 -u "$script" \
    "$root" > output.out 2> run.err || status=$?
else
  export OMP_PLACES=cores OMP_PROC_BIND=close XMVB_NEO_TRACE=scalar
  unset XMVB_NEO_CANDIDATE_AUDIT
  binary="$root/shell-screening-fix-1253607/xmvb-cpp.exe"
  steps="${tag#neo-}"
  sha256sum "$binary" input.xmi cc-pvdz.gbs > hashes.txt
  srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    /usr/bin/time -v -o resource.time "$binary" input.xmi \
    --optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal \
    --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 \
    --max-iterations "$steps" --neo-trace neo.tsv \
    --dump-final-orbital-value-table-bin final-raw.bin > output.out 2> run.err || status=$?
fi
printf '%s\n' "$status" > run.exit
# Native exit 2 is expected for the deliberately bounded path recording.
if [[ "$tag" == neo-* && "$status" == 2 ]]; then exit 0; fi
exit "$status"
