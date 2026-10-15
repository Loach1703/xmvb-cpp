#!/usr/bin/env bash
#SBATCH --job-name=fecl2-keyframe-curvature
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=fecl2-keyframe-curvature-%j.out
#SBATCH --error=fecl2-keyframe-curvature-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
set -u
root="$SLURM_SUBMIT_DIR"
run="$root/fecl2-keyframe-curvature-$SLURM_JOB_ID"
mkdir -p "$run"
cp "$root/fecl2-ritz-cost-1248947/xmvb-cpp.exe" "$run/"
cp "$root/fecl2-ritz-cost-1248947/native/input.xmi" "$run/"
cp "$root/fecl2-ritz-cost-1248947/native/def2-tzvp.gbs" "$run/"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export XMVB_NEO_TRACE=1 XMVB_NEO_CANDIDATE_AUDIT=1
cd "$run"
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
sha256sum xmvb-cpp.exe input.xmi def2-tzvp.gbs > hashes.txt
status=0
# Existing bounded audit: no full Hessian; the 419-coordinate metric is
# independently certified within the audit's fixed 16 MiB allocation limit.
# Extra keyframe evaluations are diagnostics, not a native speed measurement.
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
 ./xmvb-cpp.exe input.xmi --optimizer-backend neo --eigensolver davidson \
 --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
 --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
[[ "$status" == 0 || "$status" == 2 ]]
grep -q 'ciah_keyframe_direction_audit' run.err
