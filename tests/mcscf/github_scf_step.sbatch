#!/usr/bin/env bash
#SBATCH --job-name=github-scf-step
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --sockets-per-node=1
#SBATCH --cores-per-socket=16
#SBATCH --threads-per-core=1
#SBATCH --mem=48G
#SBATCH --time=00:20:00
#SBATCH --output=github-scf-step-%j.out
#SBATCH --error=github-scf-step-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
if [[ ${1:-} == --build-baseline && -d ${2:-} ]]; then
 cmake --build "$2" --target xmvb -j8
 exit 0
fi
[[ -x ${2:-} && -x ${3:-} ]] || exit 64
[[ $1 == COV || $1 == FULL || $1 == SINGLE ]] || exit 64
source_dir="$SLURM_SUBMIT_DIR"
space=$1
baseline=$(realpath "$2")
current=$(realpath "$3")
artifact="$source_dir/github-scf-step-$SLURM_JOB_ID"
mode=${4:-paired}
[[ $mode == paired || $mode == native || $mode == phases || $mode == vtune || $mode == converge ]] || exit 64
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=8 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
unset XMVB_PHASE_TIMING XMVB_NEO_TRACE
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
env | sort > "$artifact/environment.txt"
printf '%s\n' "$space" > "$artifact/space.txt"
printf '%s\n' '946934570c096ba5072923effa7fb063736731f1' > "$artifact/baseline-commit.txt"
printf '%s\n' "${5:-03919978 (same numerical source as e5ef2fe1)}" > "$artifact/current-commit.txt"
variants=(current)
[[ $mode == paired ]] && variants+=(baseline)
for variant in "${variants[@]}"; do
 run_dir="$artifact/$variant"
 mkdir -p "$run_dir"
 binary=$current
 if [[ $variant == baseline ]]; then binary=$baseline; fi
 cp "$binary" "$run_dir/xmvb-cpp.exe"
 # Identical original geometry, HAO supports and initial orbital coefficients.
 sed -E -e 's/^iscf=.*/iscf=9/I' -e "s/^str=.*/str=$space/I" \
  "$source_dir/testdata/vbscf/csd_LOFLEA_vb.xmi" > "$run_dir/input.xmi"
 if [[ $space == SINGLE ]]; then
  # LOFLEA has 68 electrons: 28 inactive pairs and 12 active electrons.
  sed -i '/^str=/Id' "$run_dir/input.xmi"
  printf '\n$STR\n1:28 29 30 31 32 33 34 35 36 37 38 39 40\n$END\n' >> "$run_dir/input.xmi"
 fi
 cp "$source_dir/basis/cc-pvdz.gbs" "$run_dir/"
 (
  cd "$run_dir"
  sha256sum input.xmi cc-pvdz.gbs xmvb-cpp.exe > hashes.sha256
  launcher=()
  iterations=(--max-iterations 1)
  [[ $mode == converge ]] && iterations=()
  if [[ $mode == phases || $mode == vtune ]]; then
   export XMVB_PHASE_TIMING=1 XMVB_NEO_TRACE=scalar
  fi
  if [[ $mode == vtune ]]; then
   profiler=/opt/intel/oneapi/vtune/latest/bin64/vtune
   [[ -x $profiler ]] || exit 69
   "$profiler" -version > profiler-version.txt 2>&1
   launcher=("$profiler" -collect hotspots -knob sampling-mode=sw -result-dir "$run_dir/capture.vtune" --)
  fi
  set +e
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
   /usr/bin/time -v -o resource.time \
   python3 "$source_dir/tests/mcscf/pin_ci.py" 8 "$artifact/affinity.jsonl" \
   "${launcher[@]}" "$run_dir/xmvb-cpp.exe" input.xmi --optimizer-backend neo --eigensolver davidson \
   --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 \
   --energy-tolerance 1e-7 "${iterations[@]}" --neo-trace neo.tsv > output.out 2> run.err
  status=$?
  set -e
  printf '%s\n' "$status" > run.exit
  if [[ $mode == vtune ]]; then
   for report in hotspots summary; do
    "$profiler" -report "$report" -r "$run_dir/capture.vtune" \
     > "vtune-$report.txt" 2> "vtune-$report.err"
   done
  fi
 )
done
# Application exit 2 is expected only if output confirms max_iterations;
# preserve all other failures and never report them as a completed SCF step.
