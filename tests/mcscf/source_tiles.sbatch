#!/usr/bin/env bash
#SBATCH --job-name=source-tiles
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:45:00
#SBATCH --output=source-tiles-%j.out
#SBATCH --error=source-tiles-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
deps=/home/guqqgroup/xiatao/software/xmvb-deps-20261012
artifact="$source_dir/source-tiles-$SLURM_JOB_ID"
mkdir -p "$artifact"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="$deps/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_CPUS_PER_TASK=)' > "$artifact/environment.txt"
cmake -S "$source_dir" -B "$source_dir/build" -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran \
 '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG -g1 -fno-omit-frame-pointer' \
 -DCINT_ROOT_DIR="$deps" -DLAPACK_ROOT_DIR="$deps" -DEIGEN3_ROOT_DIR="$deps" \
 -DXMVB_BUILD_DEV_TARGETS=ON -DXMVB_ENABLE_EIGEN_BLAS_BACKEND=ON \
 > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build" --target xmvb test_structure_tiles \
 test_orthogonal_direct_ci test_coupled_structure test_product_frontier \
 test_representation_overlap test_ci_transport test_selected_ci_transport \
 test_selected_metric -j8 > "$artifact/build.log" 2>&1
sha256sum "$source_dir/build/src/test_structure_tiles" \
 "$source_dir/build/src/xmvb-cpp.exe" > "$artifact/binaries.sha256"
sha256sum "$source_dir/src/mcscf/ci/spin_channels.cpp" \
 "$source_dir/src/mcscf/ci/sigma.hpp" \
 "$source_dir/src/vbscf/structures/assembly/source_tiles.cpp" \
 "$source_dir/src/vbscf/structures/assembly/sources.hpp" \
 "$source_dir/src/vbscf/structures/assembly/action.cpp" \
 "$source_dir/src/vbscf/structures/assembly/action.hpp" \
 "$source_dir/src/vbscf/structures/orthogonal_ci/exterior_transform.cpp" \
 "$source_dir/src/vbscf/structures/orthogonal_ci/exterior_transform.hpp" \
 "$source_dir/tests/vbscf/unit/test_structure_tiles.cpp" > "$artifact/sources.sha256"
srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores env OMP_NUM_THREADS=2 \
 ctest --test-dir "$source_dir/build" \
 -R '^(orthogonal_direct_ci|coupled_structure|structure_tiles|product_frontier|representation_overlap|ci_transport|selected_ci_transport|selected_metric)$' \
 --output-on-failure > "$artifact/ctest.out" 2>&1
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 bash -c 'grep Cpus_allowed_list /proc/self/status' > "$artifact/affinity.txt"
for count in 1 3 20 132; do
 for mode in selected12 selected12-parity sources12; do
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
   /usr/bin/time -v -o "$artifact/$count-$mode.time" \
   "$source_dir/build/src/test_structure_tiles" --"$mode" "$count" \
   > "$artifact/$count-$mode.out" 2> "$artifact/$count-$mode.err"
 done
done
for orbitals in 12 14; do
 for mode in source-size source-size-parity; do
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
   /usr/bin/time -v -o "$artifact/$orbitals-$mode.time" \
   "$source_dir/build/src/test_structure_tiles" --"$mode" "$orbitals" \
   > "$artifact/$orbitals-$mode.out" 2> "$artifact/$orbitals-$mode.err"
 done
done
# Keep stage timings and optional software sampling separate from native runs.
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 env XMVB_PHASE_TIMING=1 "$source_dir/build/src/test_structure_tiles" --source-size 14 \
 > "$artifact/phases.out" 2> "$artifact/phases.err"
profiler=/opt/oneapi/2025.2.1.44/vtune/2025.4/bin64/vtune
if [[ -x $profiler ]]; then
 srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
  "$profiler" -collect hotspots -knob sampling-mode=sw -result-dir "$artifact/capture" -- \
  "$source_dir/build/src/test_structure_tiles" --source-size 14 \
  > "$artifact/capture.out" 2> "$artifact/capture.err"
 srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
  "$profiler" -report hotspots -r "$artifact/capture" > "$artifact/hotspots.txt" 2> "$artifact/report.err"
else
 printf '%s\n' "VTune unavailable at $profiler on compute node" > "$artifact/profiler.txt"
fi
