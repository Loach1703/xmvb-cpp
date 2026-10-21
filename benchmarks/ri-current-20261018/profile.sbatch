#!/usr/bin/env bash
#SBATCH --job-name=RI-current-profile
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=32G
#SBATCH --time=00:10:00
#SBATCH --output=profile-%j.out
#SBATCH --error=profile-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# -ge 1 && $# -le 2 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0
set -u
frozen=$(realpath -e "$1")
profiler=${2:-vtune}
[[ $profiler == vtune || $profiler == gprofng ]] || exit 64
source_dir=$(realpath -e "$SLURM_SUBMIT_DIR/../..")
artifact="$SLURM_SUBMIT_DIR/profile-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LD_LIBRARY_PATH=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle/lib:${LD_LIBRARY_PATH:-}
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
vtune=/opt/oneapi/2025.2.1.44/vtune/2025.4/bin64/vtune
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
if [[ $profiler == vtune ]]; then
  "$vtune" -version > "$artifact/profiler-version.txt" 2>&1
else
  gprofng --version > "$artifact/profiler-version.txt" 2>&1
fi
for label in c12-cov FeCl2; do
  run_dir="$artifact/$label"
  mkdir -p "$run_dir"
  if [[ $label == c12-cov ]]; then
    cp "$source_dir/benchmarks/c12-c14-channels-20261018/cases/c12-cov/input.xmi" "$run_dir/"
  else
    cp "$source_dir/testdata/vbscf/FeCl2.xmi" "$run_dir/input.xmi"
  fi
  cp "$source_dir/basis/"{cc-pvdz.gbs,cc-pvdz-jkfit.gbs,def2-tzvp.gbs,def2-tzvp-jkfit.gbs} "$run_dir/"
  cd "$run_dir"
  sha256sum "$frozen/bin/xmvb-cpp.exe" input.xmi *.gbs > hashes.txt
  export RI_PROFILE_BINARY="$frozen/bin/xmvb-cpp.exe" RI_PROFILE_DIR="$run_dir" RI_PROFILE_TOOL="$vtune"
  export RI_PROFILE_MODE="$profiler"
  srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
cd "$RI_PROFILE_DIR"
grep Cpus_allowed_list /proc/self/status > affinity.txt
if [[ $RI_PROFILE_MODE == vtune ]]; then
 "$RI_PROFILE_TOOL" -collect hotspots -knob sampling-mode=sw -result-dir capture -- \
  "$RI_PROFILE_BINARY" input.xmi --optimizer-backend ciah --standard-two-electron-mode ri \
  --ci-workspace-mib 4096 --ciah-trace ciah.tsv > output.out 2> run.err
"$RI_PROFILE_TOOL" -report summary -r capture > summary.txt 2> summary.err
"$RI_PROFILE_TOOL" -report hotspots -r capture -group-by function -format csv > hotspots.csv 2> hotspots.err
 "$RI_PROFILE_TOOL" -report hotspots -r capture -group-by thread -format csv > threads.csv 2> threads.err
else
 gprofng collect app -p on -o capture.er "$RI_PROFILE_BINARY" input.xmi \
   --optimizer-backend ciah --standard-two-electron-mode ri \
   --ci-workspace-mib 4096 --ciah-trace ciah.tsv > output.out 2> run.err
 gprofng display text -header -statistics capture.er > quality.txt 2> quality.err
 gprofng display text -limit 80 -functions capture.er > functions.txt 2> functions.err
 gprofng display text -limit 100 -calltree capture.er > calltree.txt 2> calltree.err
fi
RUN
done
