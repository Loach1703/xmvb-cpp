#!/usr/bin/env bash
#SBATCH --job-name=RI-current-native
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=00:30:00
#SBATCH --output=ri-%j.out
#SBATCH --error=ri-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# -ge 1 && $# -le 2 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0
set -u
frozen=$(realpath -e "$1")
mode=${2:-native}
case "$mode" in native|phase|phase-fe|exact-c12|exact-quick) ;; *) exit 64 ;; esac
source_dir=$(realpath -e "$SLURM_SUBMIT_DIR/../..")
artifact="$SLURM_SUBMIT_DIR/$mode-$SLURM_JOB_ID"
mkdir -p "$artifact/bin"
cp "$frozen/bin/xmvb-cpp.exe" "$artifact/bin/"
cp "$frozen/source.sha256" "$frozen/CMakeCache.txt" "$artifact/"
cp "$SLURM_SUBMIT_DIR/run.sbatch" "$artifact/"
export LD_LIBRARY_PATH=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle/lib:${LD_LIBRARY_PATH:-}
export LC_ALL=C OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export XMVB_CIAH_TRACE=scalar
unset XMVB_PHASE_TIMING
if [[ $mode == phase || $mode == phase-fe ]]; then export XMVB_PHASE_TIMING=1; fi
labels=(241 c12-cov c14-cov FeCl2 MnF2)
if [[ $mode == phase ]]; then labels=(c12-cov); fi
if [[ $mode == phase-fe ]]; then labels=(FeCl2); fi
integrals=ri
if [[ $mode == exact-c12 ]]; then labels=(c12-cov); integrals=exact; fi
if [[ $mode == exact-quick ]]; then labels=(241 FeCl2 MnF2); integrals=exact; fi
export RI_BENCH_INTEGRALS="$integrals"
for label in "${labels[@]}"; do
  case_dir="$artifact/$label"
  mkdir -p "$case_dir"
  case "$label" in
    c12-cov|c14-cov)
      input="$source_dir/benchmarks/c12-c14-channels-20261018/cases/$label/input.xmi" ;;
    241) input="$source_dir/testdata/vbscf/241_vb.xmi" ;;
    *) input="$source_dir/testdata/vbscf/$label.xmi" ;;
  esac
  cp "$input" "$case_dir/input.xmi"
  cp "$source_dir/basis/"{cc-pvdz.gbs,cc-pvdz-jkfit.gbs,def2-tzvp.gbs,def2-tzvp-jkfit.gbs} "$case_dir/"
  cd "$case_dir"
  hostname > node.txt
  lscpu > cpu.txt
  module list > modules.txt 2>&1
  env | sort > environment.txt
  sha256sum "$artifact/bin/xmvb-cpp.exe" input.xmi *.gbs > hashes.txt
  printf '%s\n' "$SLURM_JOB_ID" > slurm-job.txt
  export RI_BENCH_DIR="$case_dir" RI_BENCH_BINARY="$artifact/bin/xmvb-cpp.exe"
  # Each case is a separate native process, sharing only the reserved core mask.
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores bash -se <<'RUN'
cd "$RI_BENCH_DIR"
grep Cpus_allowed_list /proc/self/status > affinity.txt
set +e
/usr/bin/time -v -o resource.time "$RI_BENCH_BINARY" input.xmi \
  --optimizer-backend ciah --standard-two-electron-mode "$RI_BENCH_INTEGRALS" \
  --ci-workspace-mib 4096 --ciah-trace ciah.tsv > output.out 2> run.err
status=$?
printf '%s\n' "$status" > run.exit
tail -n90 output.out > report-tail.txt
exit "$status"
RUN
done
