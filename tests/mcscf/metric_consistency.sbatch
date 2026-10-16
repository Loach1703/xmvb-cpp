#!/usr/bin/env bash
#SBATCH --job-name=metric-consistency
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:20:00
#SBATCH --output=metric-consistency-%j.out
#SBATCH --error=metric-consistency-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -d ${1:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
point_dir=$(realpath "$1")
artifact="$source_dir/metric-consistency-$SLURM_JOB_ID"
mkdir -p "$artifact/point-input"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_CPUS_PER_TASK=)' > "$artifact/environment.txt"
cmake -S "$source_dir" -B "$source_dir/build" > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build" --target test_metric_point test_sparse_orbital_quotient \
 test_ciah test_shifted_metric_preconditioner test_structure_tiles test_coupled_structure \
 test_orthogonal_direct_ci test_support_preserving_gauge test_keyframe -j8 > "$artifact/build.log" 2>&1
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
sha256sum "$source_dir/src/vbscf/orbitals/charts/physical_metric."{cpp,hpp} \
 "$source_dir/src/vbscf/optimization/backends/ciah.cpp" \
 "$source_dir/src/vbscf/derivatives/gradient/active_space/evaluator.cpp" \
 "$source_dir/tests/vbscf/diagnostics/metric_point.cpp" > "$artifact/sources.sha256"
sha256sum "$source_dir/build/src/test_metric_point" > "$artifact/binary.sha256"
srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores env OMP_NUM_THREADS=2 \
 ctest --test-dir "$source_dir/build" --output-on-failure \
 -R '^(sparse_orbital_quotient|ciah|shifted_metric_preconditioner|structure_tiles|coupled_structure|orthogonal_direct_ci|support_preserving_gauge|keyframe)$' \
 > "$artifact/ctest.out" 2>&1
cp "$point_dir/input.xmi" "$point_dir/cc-pvdz.gbs" "$point_dir/orbitals.bin" "$artifact/point-input/"
cd "$artifact/point-input"
sha256sum input.xmi cc-pvdz.gbs orbitals.bin > input.sha256
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 /usr/bin/time -v -o "$artifact/audit.time" "$source_dir/build/src/test_metric_point" \
 input.xmi orbitals.bin "$artifact/point" > "$artifact/audit.out" 2> "$artifact/audit.err"
cd "$source_dir"
# Same native capture, with no scalar events; this is a bounded diagnostic,
# not a claim of convergence at the diagnostic iteration limit.
if [[ ${2:-} == --span-regression ]]; then
 cmake --build "$source_dir/build" --target xmvb -j8 > "$artifact/application-build.log" 2>&1
 sha256sum "$source_dir/build/src/xmvb-cpp.exe" >> "$artifact/binary.sha256"
 mkdir -p "$artifact/FeCl2"
 cp "$source_dir/testdata/vbscf/FeCl2.xmi" "$artifact/FeCl2/input.xmi"
 cp "$source_dir/basis/def2-tzvp.gbs" "$artifact/FeCl2/"
 cd "$artifact/FeCl2"
 sha256sum input.xmi def2-tzvp.gbs > input.sha256
 srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
  /usr/bin/time -v -o resource.time "$source_dir/build/src/xmvb-cpp.exe" input.xmi \
  --optimizer-backend ciah --eigensolver davidson --orbital-preconditioner hessian-diagonal \
  --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --ciah-trace neo.tsv \
  > output.out 2> run.err
elif [[ ${2:-} == --regression ]]; then
 bash "$source_dir/tests/mcscf/source_cover.sbatch" --production "$3"
else
 bash "$source_dir/tests/mcscf/selected_metric.sbatch" "$point_dir"
fi
