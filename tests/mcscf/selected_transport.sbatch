#!/usr/bin/env bash
#SBATCH --job-name=selected-transport
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:20:00
#SBATCH --output=selected-transport-%j.out
#SBATCH --error=selected-transport-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
artifact="$source_dir/selected-transport-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_CPUS_PER_TASK=)' > "$artifact/environment.txt"
sha256sum "$source_dir/src/vbscf/optimization/coupled/ci_transport."{cpp,hpp} \
 "$source_dir/tests/vbscf/unit/test_selected_metric.cpp" > "$artifact/sources.sha256"
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 bash -c 'grep Cpus_allowed_list /proc/self/status' > "$artifact/affinity.txt"
if [[ ${1:-} == --build ]]; then
 cmake -S "$source_dir" -B "$source_dir/build" > "$artifact/configure.log" 2>&1
 cmake --build "$source_dir/build" --target xmvb test_ci_factors test_structure_tiles \
  test_orthogonal_direct_ci test_coupled_structure test_product_frontier \
  test_representation_overlap test_ci_transport test_selected_ci_transport \
  test_selected_metric test_eigen_response test_keyframe test_davidson -j8 > "$artifact/build.log" 2>&1
 cp "$source_dir/build/CMakeCache.txt" "$artifact/"
 sha256sum "$source_dir/build/src/xmvb-cpp.exe" > "$artifact/binary.sha256"
 srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores env OMP_NUM_THREADS=2 \
  ctest --test-dir "$source_dir/build" --output-on-failure \
  -R '^(ci_factors|orthogonal_direct_ci|coupled_structure|structure_tiles|product_frontier|representation_overlap|ci_transport|selected_ci_transport|selected_metric|eigen_response|keyframe|davidson)$' \
  > "$artifact/ctest.out" 2>&1
 cp "$source_dir/build/Testing/Temporary/LastTest.log" "$artifact/test-details.log"
 exit 0
fi
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/"
sha256sum "$artifact/xmvb-cpp.exe" > "$artifact/binary.sha256"
run_case() {
 local case_dir=$1
 cd "$case_dir"
 sha256sum input.xmi *.gbs > input.sha256
 set +e
 srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
  /usr/bin/time -v -o resource.time "$artifact/xmvb-cpp.exe" input.xmi \
  --optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal \
  --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --neo-trace neo.tsv > output.out 2> run.err
 local status=$?
 set -e
 printf '%s\n' "$status" > run.exit
 return "$status"
}
if [[ ${1:-} == --240 && -d ${2:-} ]]; then
 bundle=$(realpath "$2")
 labels=(cov ion01 ion02 ion03 full)
 keywords=('COV' 'ION(0:1)' 'ION(0:2)' 'ION(0:3)' 'FULL')
 for index in 0 1 2 3 4; do
  mkdir -p "$artifact/native/${labels[index]}"
  sed -E -e "s/^str=.*/str=${keywords[index]}/I" -e 's/^iscf=.*/iscf=9/I' \
   "$bundle/240_vb.xmi" > "$artifact/native/${labels[index]}/input.xmi"
  cp "$bundle/cc-pvdz.gbs" "$artifact/native/${labels[index]}/"
  run_case "$artifact/native/${labels[index]}"
 done
elif [[ ${1:-} == --FeCl2 ]]; then
 mkdir -p "$artifact/FeCl2"
 cp "$source_dir/testdata/vbscf/FeCl2.xmi" "$artifact/FeCl2/input.xmi"
 cp "$source_dir/basis/def2-tzvp.gbs" "$artifact/FeCl2/"
 run_case "$artifact/FeCl2"
else
 exit 64
fi
