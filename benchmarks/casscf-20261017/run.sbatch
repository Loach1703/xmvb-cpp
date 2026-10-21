#!/usr/bin/env bash
#SBATCH --job-name=casscf-match
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --output=casscf-%j.out
#SBATCH --error=casscf-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# -eq 2 ]] || exit 64
source_dir=$(realpath -e "$1")
case_name=$2
[[ $case_name == F2 || $case_name == FeCl2 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0 anaconda3/25.06
set -u
envprefix=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
bundle=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle
wheel=/home/guqqgroup/taoxia/xmvb-runs/oeo-pyscf-20261011/wheels/pyscf-2.14.0-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
python_bin=/opt/anaconda3/25.06/bin/python3
binary="$source_dir/build/src/casscf_fixture"
fecl2_layout=/home/guqqgroup/taoxia/xmvb-runs/oeo-pyscf-20261011/runs-282255/FeCl2/cpp/input.molden
artifact="$SLURM_SUBMIT_DIR/casscf-${case_name}-${SLURM_JOB_ID}"
mkdir -p "$artifact"
[[ -f "$source_dir/parent-revision.txt" && -f "$source_dir/source.sha256" && -x "$binary" ]] || {
  echo "qualified source manifest or shared casscf_fixture binary is missing" >&2
  exit 66
}
export LD_LIBRARY_PATH="$bundle/lib:$envprefix/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export PYTHONUNBUFFERED=1
hostname > "$artifact/node.txt"
lscpu > "$artifact/lscpu.txt"
module list > "$artifact/modules.txt" 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > "$artifact/environment.txt"
cp "$source_dir/parent-revision.txt" "$artifact/"
cp "$source_dir/source.sha256" "$artifact/"
(cd "$source_dir" && sha256sum --check source.sha256) \
  > "$artifact/source-verify.out" 2> "$artifact/source-verify.err"
"$python_bin" --version > "$artifact/python-version.txt" 2>&1
g++ --version > "$artifact/compiler-version.txt"
mkdir "$artifact/python-overlay"
"$python_bin" -m zipfile -e "$wheel" "$artifact/python-overlay"
export PYTHONPATH="$artifact/python-overlay"
"$python_bin" -c 'import numpy, scipy, h5py, pyscf; print(numpy.__version__, scipy.__version__, h5py.__version__, pyscf.__version__)' \
  > "$artifact/pyscf-version.txt"
sha256sum "$binary" "$source_dir/benchmarks/casscf-20261017/"*.py \
  "$source_dir/basis/cc-pvdz.gbs" "$source_dir/basis/def2-tzvp.gbs" "$wheel" \
  > "$artifact/input-and-binary.sha256"
if [[ $case_name == FeCl2 ]]; then
  [[ -r "$fecl2_layout" ]] || { echo "qualified FeCl2 Molden layout is missing" >&2; exit 66; }
  sha256sum "$fecl2_layout" >> "$artifact/input-and-binary.sha256"
fi
ldd "$binary" > "$artifact/dependencies.txt"
! grep -q 'not found' "$artifact/dependencies.txt"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  bash -c 'grep Cpus_allowed_list /proc/self/status' > "$artifact/affinity.txt"

scripts="$source_dir/benchmarks/casscf-20261017"
export_options=()
if [[ $case_name == FeCl2 ]]; then
  export_options=(--fecl2-molden-layout "$fecl2_layout")
fi
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o "$artifact/export.time" \
  env -u OMP_PLACES OMP_PROC_BIND=false "$python_bin" -u "$scripts/export_fixture.py" \
  "$case_name" "$source_dir" "$artifact/fixture" "${export_options[@]}" \
  > "$artifact/export.out" 2> "$artifact/export.err"
sha256sum "$artifact/fixture/ao-mo.bin" "$artifact/fixture/metadata.json" \
  >> "$artifact/input-and-binary.sha256"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o "$artifact/initial-cpp.time" \
  "$binary" --initial-only "$artifact/fixture/ao-mo.bin" \
  > "$artifact/initial-cpp.out" 2> "$artifact/initial-cpp.err"
"$python_bin" -u "$scripts/check_initial.py" \
  "$artifact/fixture/metadata.json" "$artifact/initial-cpp.out" \
  > "$artifact/initial-parity.out"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o "$artifact/pyscf.time" \
  env -u OMP_PLACES OMP_PROC_BIND=false "$python_bin" -u "$scripts/pyscf_newton.py" \
  "$artifact/fixture/ao-mo.bin" "$artifact/fixture/metadata.json" \
  "$source_dir" "$artifact/pyscf" \
  > "$artifact/pyscf.out" 2> "$artifact/pyscf.err"
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  env XMVB_CIAH_TRACE=scalar /usr/bin/time -v -o "$artifact/cpp.time" \
  "$binary" "$artifact/fixture/ao-mo.bin" "$artifact/cpp-trajectory.csv" \
  > "$artifact/cpp.out" 2> "$artifact/cpp.err" || cpp_status=$?
printf '%s\n' "${cpp_status:-0}" > "$artifact/cpp.exit"
if [[ ${cpp_status:-0} != 0 ]]; then
  echo "C++ CASSCF failed; PySCF result and C++ diagnostic logs are retained" >&2
  exit "$cpp_status"
fi
"$python_bin" -u "$scripts/compare.py" \
  "$artifact/fixture/metadata.json" "$artifact/pyscf/result.json" \
  "$artifact/cpp.out" "$artifact/comparison.json" > "$artifact/compare.out"
