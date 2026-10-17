#!/usr/bin/env bash
#SBATCH --job-name=C12-spaces
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=64G
#SBATCH --output=C12-spaces-%A_%a.out
#SBATCH --error=C12-spaces-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
set -u
source_dir="$SLURM_SUBMIT_DIR"
export LC_ALL=C OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE

if [[ ${1:-} == --prepare && $# == 3 ]]; then
  qualified=$(realpath -e "$2")
  bundle="$source_dir/c12-spaces-bundle-$SLURM_JOB_ID"
  mkdir -p "$bundle"
  grep -q '^Termination reason *: ciah_projected_gradient_tolerance' "$qualified/C14.out"
  (cd "$qualified" && sha256sum --check hashes.sha256) > "$bundle/qualification-check.txt"
  sha256sum --check "$qualified/sources.sha256" > "$bundle/source-check.txt"
  cp "$qualified/"{xmvb-cpp.exe,cc-pvdz.gbs,CMakeCache.txt,sources.sha256,ctest.out} "$bundle/"
  printf '%s\n' "$3" > "$bundle/revision.txt"
  printf '%s\n' "$qualified" > "$bundle/qualification.txt"
  cp "$source_dir/testdata/vbscf/csd_LOFLEA_vb.xmi" "$bundle/original.xmi"
  cp "$source_dir/tests/mcscf/"{c12_spaces.py,c12_spaces.sbatch,summarize_scf_step.py} "$bundle/"
  python3 "$bundle/c12_spaces.py" prepare "$bundle/original.xmi" "$bundle" > "$bundle/preparation.json"
  (cd "$bundle" && sha256sum xmvb-cpp.exe original.xmi cc-pvdz.gbs \
      c12_spaces.py c12_spaces.sbatch summarize_scf_step.py cases.json > hashes.sha256)
  exit 0
fi

[[ $# == 1 && -n ${SLURM_ARRAY_TASK_ID:-} && $SLURM_CPUS_PER_TASK == 32 ]] || exit 64
bundle=$(realpath -e "$1")
(cd "$bundle" && sha256sum --check hashes.sha256) > /dev/null
label=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["cases"][int(sys.argv[2])]["label"])' \
  "$bundle/cases.json" "$SLURM_ARRAY_TASK_ID")
artifact="$source_dir/c12-spaces-$SLURM_ARRAY_JOB_ID/$label"
mkdir -p "$artifact"
cp "$bundle/cases/$label/"{input.xmi,cc-pvdz.gbs} "$artifact/"
cp "$bundle/"{cases.json,revision.txt,qualification.txt} "$artifact/"
cd "$artifact"
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
printf '%s\n' "$bundle" > bundle.txt
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > environment.txt
sha256sum "$bundle/xmvb-cpp.exe" input.xmi cc-pvdz.gbs > hashes.sha256
set +e
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -c '
  grep Cpus_allowed_list /proc/self/status > affinity.txt
  exec /usr/bin/time -v -o resource.time "$1" input.xmi --ciah-trace ciah.tsv
' bash "$bundle/xmvb-cpp.exe" > output.out 2> run.err
status=$?
printf '%s\n' "$status" > run.exit
set -e
[[ "$status" == 0 ]]
grep -q '^Termination reason *: ciah_projected_gradient_tolerance' output.out
