#!/usr/bin/env bash
#SBATCH --job-name=C12-ionic-profile
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=64G
#SBATCH --output=C12-ionic-profile-%j.out
#SBATCH --error=C12-ionic-profile-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 2 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
set -u
bundle=$(realpath -e "$1")
label=$2
[[ $label == ion01 || $label == ion02 ]] || exit 64
artifact="$SLURM_SUBMIT_DIR/c12-ionic-profile-$SLURM_JOB_ID"
mkdir -p "$artifact"
(cd "$bundle" && sha256sum --check hashes.sha256) > "$artifact/bundle-check.txt"
cp "$bundle/cases/$label/"{input.xmi,cc-pvdz.gbs} "$artifact/"
cp "$bundle/"{revision.txt,CMakeCache.txt} "$artifact/"
# Attribution only: limit the unchanged current executable to its first macrostep.
# This capture is never a native timing or complete-convergence comparison.
python3 - "$artifact/input.xmi" <<'PY'
import pathlib, re, sys
path = pathlib.Path(sys.argv[1])
text, count = re.subn(r'^itmax=.*$', 'itmax=1', path.read_text(), flags=re.M | re.I)
assert count == 1
path.write_text(text)
PY
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export PATH="/home/guqqgroup/xiatao/software/gprofng-2.42/bin:$PATH"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/gprofng-2.42/lib:/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export XMVB_PHASE_TIMING=1 XMVB_CIAH_TRACE=scalar
export C12_PROFILE_ROOT="$artifact" C12_PROFILE_BINARY="$bundle/xmvb-cpp.exe"
hostname > "$artifact/node.txt"
module list > "$artifact/modules.txt" 2>&1
cd "$artifact"
sha256sum "$C12_PROFILE_BINARY" input.xmi cc-pvdz.gbs > hashes.sha256
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'RUN'
cd "$C12_PROFILE_ROOT"
grep Cpus_allowed_list /proc/self/status > affinity.txt
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_|XMVB_)' > environment.txt
gprofng --version > profiler-version.txt
set +e
/usr/bin/time -v -o collection.time gprofng collect app -p on -o capture.er \
    "$C12_PROFILE_BINARY" input.xmi --ciah-trace ciah.tsv > output.out 2> run.err
status=$?
printf '%s\n' "$status" > capture.exit
set -e
gprofng display text -header -statistics capture.er > quality.txt
gprofng display text -limit 60 -functions -calltree capture.er > hotspots.txt
gprofng display text -limit 60 -lines capture.er > lines.txt 2> lines.err
RUN
