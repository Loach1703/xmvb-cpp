#!/usr/bin/env bash
#SBATCH --job-name=ciah-affine
#SBATCH --partition=CPU-256C256GB
#SBATCH --qos=qos_cpu-256c256gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=affine-%j.out
#SBATCH --error=affine-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0
deps=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
runtime=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle/lib
export PATH="$deps/bin:$PATH"
export LD_LIBRARY_PATH="$runtime:$deps/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
artifact="$SLURM_SUBMIT_DIR/affine-$SLURM_JOB_ID"
mkdir -p "$artifact"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
find src tests/vbscf/unit tests/linalg/unit -type f \
  \( -name '*.cpp' -o -name '*.hpp' -o -name '*.cmake' \) -print0 |
  sort -z | xargs -0 sha256sum > "$artifact/sources.sha256"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DXMVB_BUILD_DEV_TARGETS=ON \
  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran \
  -DLAPACK_ROOT_DIR="$deps" -DEIGEN3_ROOT_DIR="$deps" -DCINT_ROOT_DIR="$deps" \
  -DLIBXC_ROOT_DIR="$deps" -DLAPACK_LIB="$runtime/libopenblas.so.0" \
  -DCINT_LIB="$runtime/libcint.so.6" > "$artifact/configure.log" 2>&1
cmake --build build --target test_coupled_structure test_eigen_response \
  test_selected_ci_transport test_ci_transport test_keyframe test_ciah xmvb \
  -j "$SLURM_CPUS_PER_TASK" > "$artifact/build.log" 2>&1
sha256sum build/src/test_coupled_structure build/src/xmvb-cpp.exe > "$artifact/binaries.sha256"
status=0
for test in coupled_structure eigen_response selected_ci_transport ci_transport keyframe ciah; do
  if srun --ntasks=1 --cpus-per-task=2 --cpu-bind=cores \
    env OMP_NUM_THREADS=1 "build/src/test_$test" > "$artifact/$test.out" 2>&1; then
    printf '%s PASS\n' "$test" >> "$artifact/results.txt"
  else
    printf '%s FAIL\n' "$test" >> "$artifact/results.txt"
    status=1
  fi
done
if [[ $status == 0 ]]; then
  for molecule in F2_OEO_PYSCF_LOW F2_SA2; do
    if srun --ntasks=1 --cpus-per-task=4 --cpu-bind=cores env OMP_NUM_THREADS=4 \
      /usr/bin/time -v -o "$artifact/$molecule.time" \
      build/src/xmvb-cpp.exe "testdata/vbscf/$molecule.xmi" \
      --optimizer-backend ciah --eigensolver davidson \
      --ciah-trace "$artifact/$molecule.tsv" > "$artifact/$molecule.out" 2>&1 &&
      grep -q 'VBSCF converged in' "$artifact/$molecule.out"; then
      printf '%s PASS\n' "$molecule" >> "$artifact/results.txt"
    else
      printf '%s FAIL\n' "$molecule" >> "$artifact/results.txt"
      status=1
    fi
  done
fi
exit "$status"
