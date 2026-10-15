#!/usr/bin/env bash
#SBATCH --job-name=casscf-common
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=00:30:00
#SBATCH --output=validate-%j.out
#SBATCH --error=validate-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5 anaconda/3/22.05
source_dir="$SLURM_SUBMIT_DIR"
deps=/home/guqqgroup/xiatao/software/xmvb-deps-20261012
artifact="$source_dir/validation-$SLURM_JOB_ID"
mkdir -p "$artifact"
export OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_MAX_ACTIVE_LEVELS=1
export LD_LIBRARY_PATH="$deps/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH=/home/guqqgroup/xiatao/xmvb-runs/shared-response-20261012/porphine-reference-audit-1247977/python
hostname > "$artifact/node.txt"
module list > "$artifact/modules.txt" 2>&1
cmake -S "$source_dir" -B "$source_dir/build" \
 -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
 -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_CXX_FLAGS_RELEASE='-O3 -DNDEBUG -g1 -fno-omit-frame-pointer' \
 -DCMAKE_PREFIX_PATH="$deps" -DLAPACK_ROOT_DIR="$deps" -DCINT_ROOT_DIR="$deps" \
 -DEIGEN3_ROOT_DIR="$deps" -DXMVB_BUILD_DEV_TARGETS=ON -DXMVB_ENABLE_CCACHE=ON \
 > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build" --target test_casci -j16 > "$artifact/casci-build.log" 2>&1
srun --exact --cpu-bind=cores "$source_dir/build/src/test_casci" > "$artifact/casci.out" 2> "$artifact/casci.err"
cd "$artifact"
srun --exact --cpu-bind=cores python3 "$source_dir/tests/mcscf/reference.py" > reference.out 2> reference.err
for name in F2 O2 H4; do
 srun --exact --cpu-bind=cores /usr/bin/time -v -o "$name.time" \
  "$source_dir/build/src/test_casci" "$name.integrals" > "$name.out" 2> "$name.err"
done
cmake --build "$source_dir/build" --target xmvb test_orthogonal_direct_ci test_davidson \
 test_ciah test_neo test_ci_transport test_representation_overlap test_selected_ci_transport \
 test_coupled_structure test_spinor_strings test_spinor_projection -j16 > vb-build.log 2>&1
srun --exact --cpu-bind=cores ctest --test-dir "$source_dir/build" \
 -R '^(casci|orthogonal_direct_ci|davidson|ciah|neo|ci_transport|representation_overlap|selected_ci_transport|coupled_structure|spinor_strings|spinor_projection|neo_f2|neo_pyscf_reference_f2|eigensolver_equivalence_f2|eigensolver_equivalence_f2_oeo)$' \
 --output-on-failure > ctest.out 2>&1
sha256sum "$source_dir/build/src/xmvb-cpp.exe" "$source_dir/build/src/test_casci" > binaries.sha256
