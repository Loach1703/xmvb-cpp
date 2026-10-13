#!/usr/bin/env bash
#SBATCH --job-name=neo-metric-reuse
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-15
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=15-00:00:00
#SBATCH --output=metric-reuse-%A_%a.out
#SBATCH --error=metric-reuse-%A_%a.err

set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
frozen="$root/fecl2-ritz-cost-${XMVB_METRIC_BUILD_JOB:?frozen build job is required}"
binary="$frozen/xmvb-cpp.exe"
index=$SLURM_ARRAY_TASK_ID
[[ "$index" != 15 ]] || index=9
row=$(sed -n "$((index + 2))p" "$root/preconditioner-suite.tsv")
IFS=$'\t' read -r case_name input_source orbital_type <<< "$row"
[[ -n "$case_name" && -f "$root/$input_source" && -x "$binary" ]] || exit 64
label="$case_name"
[[ "$SLURM_ARRAY_TASK_ID" != 15 ]] || label="$case_name-audit"
artifact="$root/metric-reuse-$SLURM_ARRAY_JOB_ID/$label"
mkdir -p "$artifact"
cp "$root/$input_source" "$artifact/input.xmi"
cp "$root/regression-1246756/241/"*.gbs "$artifact/"
cp "$root/regression-1246756/FeCl2/"*.gbs "$artifact/"
cp "$frozen/flags.txt" "$artifact/"
cd "$artifact"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
unset XMVB_NEO_TRACE XMVB_NEO_CANDIDATE_AUDIT
# Native jobs have no diagnostic-only norm evaluations. One NEW-binary capture
# attributes metric reuse separately; its wall time is not a native speed claim.
[[ "$SLURM_ARRAY_TASK_ID" != 15 ]] || export XMVB_NEO_TRACE=scalar
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_|XMVB_NEO_)' > execution-env.txt
printf '%s\n' "$orbital_type" > orbital-type.txt
printf '%s\n' 'Certified fixed-point metric inverse-image reuse; unchanged CIAH/HVP/thresholds' > source.txt
sha256sum "$binary" input.xmi *.gbs > hashes.txt
printf '%s\n' "$binary input.xmi --optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 --neo-trace neo.tsv" > command.txt
status=0
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
 "$binary" input.xmi --optimizer-backend neo --eigensolver davidson \
 --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
 --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
exit "$status"
