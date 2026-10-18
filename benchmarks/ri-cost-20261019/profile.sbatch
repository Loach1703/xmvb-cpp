#!/usr/bin/env bash
#SBATCH --job-name=RI-current-profile
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=00:10:00
#SBATCH --output=profile-%j.out
#SBATCH --error=profile-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 2 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0
set -u
frozen=$(realpath -e "$1")
label=$2
source_dir=$(realpath -e "$SLURM_SUBMIT_DIR/../..")
artifact="$SLURM_SUBMIT_DIR/profile-$SLURM_JOB_ID"
mkdir -p "$artifact/bin"
cp "$frozen/bin/xmvb-cpp.exe" "$artifact/bin/"
cp "$frozen/source.sha256" "$frozen/CMakeCache.txt" "$artifact/"
case "$label" in
  c12-cov) input="$source_dir/benchmarks/c12-c14-channels-20261018/cases/c12-cov/input.xmi" ;;
  FeCl2) input="$source_dir/testdata/vbscf/FeCl2.xmi" ;;
  *) exit 64 ;;
esac
cp "$input" "$artifact/input.xmi"
cp "$source_dir/basis/"{cc-pvdz.gbs,cc-pvdz-jkfit.gbs,def2-tzvp.gbs,def2-tzvp-jkfit.gbs} "$artifact/"
export LD_LIBRARY_PATH=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle/lib:${LD_LIBRARY_PATH:-}
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export XMVB_CIAH_TRACE=scalar XMVB_PHASE_TIMING=1 RI_PROFILE_DIR="$artifact"
cd "$artifact"
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
env | sort > environment.txt
sha256sum bin/xmvb-cpp.exe input.xmi *.gbs > hashes.txt
gprofng --version > profiler-version.txt
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
cd "$RI_PROFILE_DIR"
grep Cpus_allowed_list /proc/self/status > affinity.txt
set +e
/usr/bin/time -v -o resource.time gprofng collect app -p on -o capture.er \
 ./bin/xmvb-cpp.exe input.xmi --optimizer-backend ciah \
 --standard-two-electron-mode ri --ci-workspace-mib 4096 --ciah-trace ciah.tsv \
 > output.out 2> run.err
status=$?
printf '%s\n' "$status" > run.exit
gprofng display text -header -statistics capture.er > quality.txt 2> quality.err
gprofng display text -limit 70 -functions capture.er > functions.txt 2> functions.err
gprofng display text -limit 140 -calltree capture.er > calltree.txt 2> calltree.err
gprofng display text -limit 40 -threads capture.er > threads.txt 2> threads.err
exit "$status"
RUN
