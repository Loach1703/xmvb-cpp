#!/usr/bin/env bash
#SBATCH --job-name=remaining-cost
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=00:15:00
#SBATCH --output=remaining-cost-%A_%a.out
#SBATCH --error=remaining-cost-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
frozen="$root/fecl2-ritz-cost-1250646"
binary="$frozen/xmvb-cpp.exe"
variant=scalar
[[ "$SLURM_ARRAY_TASK_ID" == 0 ]] || variant="${XMVB_COST_CAPTURE:-sampled}"
artifact="$root/remaining-cost-$SLURM_ARRAY_JOB_ID/LOFLEA-$variant"
mkdir -p "$artifact"
cp "$root/keyframe-lean-1250648/LOFLEA/input.xmi" "$artifact/"
cp "$root/keyframe-lean-1250648/LOFLEA/"*.gbs "$artifact/"
cp "$frozen/flags.txt" "$artifact/"
cd "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_NEO_TRACE XMVB_NEO_CANDIDATE_AUDIT
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > execution-env.txt
printf '%s\n' 'b9d775e4: unchanged new binary frozen in build 1250646; no old binary rerun' > source.txt
args=(input.xmi --optimizer-backend neo --eigensolver davidson
 --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3
 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 --neo-trace neo.tsv)
if [[ "$variant" == scalar ]]; then
 export XMVB_NEO_TRACE=scalar
else
 sed -i -E 's/itmax=[0-9]+/itmax=2/g' input.xmi
 if [[ "$variant" == phases ]]; then
  export XMVB_PHASE_TIMING=1 XMVB_NEO_TRACE=scalar
 else
  vtune=/opt/intel/oneapi/vtune/latest/bin64/vtune
  "$vtune" -version > profiler-version.txt 2>&1
 fi
fi
sha256sum "$binary" input.xmi *.gbs > hashes.txt
printf '%s ' "$binary" "${args[@]}" > command.txt
status=0
if [[ "$variant" != sampled ]]; then
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
  "$binary" "${args[@]}" > output.out 2> run.err || status=$?
else
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  "$vtune" -collect hotspots -knob sampling-mode=sw -result-dir capture \
  -- "$binary" "${args[@]}" > output.out 2> run.err || status=$?
fi
printf '%s\n' "$status" > run.exit
if [[ "$variant" == sampled && -d capture ]]; then
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores "$vtune" -report hotspots -r capture \
  -format csv -csv-delimiter comma > hotspots.csv 2> report.err
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores "$vtune" -report callstacks -r capture \
  -format csv -csv-delimiter comma > callstacks.csv 2>> report.err
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores "$vtune" -report summary -r capture \
  > summary.txt 2>> report.err
fi
[[ "$status" == 0 || ( "$variant" != scalar && "$status" == 2 ) ]]
