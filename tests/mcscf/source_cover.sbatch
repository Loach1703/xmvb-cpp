#!/usr/bin/env bash
#SBATCH --job-name=source-cover
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:20:00
#SBATCH --output=source-cover-%j.out
#SBATCH --error=source-cover-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_dir="$SLURM_SUBMIT_DIR"
deps=/home/guqqgroup/xiatao/software/xmvb-deps-20261012
artifact="$source_dir/source-cover-$SLURM_JOB_ID"
mkdir -p "$artifact"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="$deps/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_CPUS_PER_TASK=)' > "$artifact/environment.txt"
cmake -S "$source_dir" -B "$source_dir/build" > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build" --target test_ci_factors test_structure_tiles test_orthogonal_direct_ci \
 test_coupled_structure test_product_frontier test_representation_overlap test_ci_transport \
 test_selected_ci_transport test_restricted_metric -j8 > "$artifact/build.log" 2>&1
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
sha256sum "$source_dir/build/src/test_structure_tiles" > "$artifact/binary.sha256"
sha256sum "$source_dir/build/src/test_ci_factors" >> "$artifact/binary.sha256"
sha256sum "$source_dir/src/vbscf/structures/assembly/"{source_tiles.cpp,sources.cpp,sources.hpp,action.cpp} \
 "$source_dir/tests/vbscf/unit/test_structure_tiles.cpp" > "$artifact/sources.sha256"
sha256sum "$source_dir/src/mcscf/ci/"{factors.cpp,factors.hpp,factors_internal.hpp,factor_generator.cpp,sigma.cpp,sigma.hpp,spin_channels.cpp} >> "$artifact/sources.sha256"
srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores env OMP_NUM_THREADS=2 \
 ctest --test-dir "$source_dir/build" \
 -R '^(ci_factors|orthogonal_direct_ci|coupled_structure|structure_tiles|product_frontier|representation_overlap|ci_transport|selected_ci_transport|restricted_metric)$' \
 --output-on-failure > "$artifact/ctest.out" 2>&1
[[ ${1:-} == --tests-only ]] && exit 0
if [[ ${1:-} == --generator16 ]]; then
 srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
  /usr/bin/time -v -o "$artifact/16-factor-generator.time" \
  "$source_dir/build/src/test_structure_tiles" --factor-generator 16 \
  > "$artifact/16-factor-generator.out" 2> "$artifact/16-factor-generator.err"
 exit 0
fi
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 bash -c 'grep Cpus_allowed_list /proc/self/status' > "$artifact/affinity.txt"
for n in 12 14; do
 for mode in source-size source-cover-size source-cover-reference; do
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
   /usr/bin/time -v -o "$artifact/$n-$mode.time" \
   "$source_dir/build/src/test_structure_tiles" --"$mode" "$n" \
   > "$artifact/$n-$mode.out" 2> "$artifact/$n-$mode.err"
 done
done
for n in 12 14; do
 for mode in factor-adjoint factor-adjoint-reference factor-direction factor-direction-reference factor-generator factor-generator-reference; do
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
   /usr/bin/time -v -o "$artifact/$n-$mode.time" \
   "$source_dir/build/src/test_structure_tiles" --"$mode" "$n" \
   > "$artifact/$n-$mode.out" 2> "$artifact/$n-$mode.err"
 done
done
for count in 1 3 20 132; do
 srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
  /usr/bin/time -v -o "$artifact/$count-selected12.time" \
  "$source_dir/build/src/test_structure_tiles" --selected12 "$count" \
  > "$artifact/$count-selected12.out" 2> "$artifact/$count-selected12.err"
done
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 env XMVB_PHASE_TIMING=1 "$source_dir/build/src/test_structure_tiles" --source-size 14 \
 > "$artifact/phases.out" 2> "$artifact/phases.err"
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 env XMVB_PHASE_TIMING=1 "$source_dir/build/src/test_structure_tiles" --factor-adjoint 14 \
 > "$artifact/adjoint-phases.out" 2> "$artifact/adjoint-phases.err"
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 env XMVB_PHASE_TIMING=1 "$source_dir/build/src/test_structure_tiles" --factor-direction 14 \
 > "$artifact/direction-phases.out" 2> "$artifact/direction-phases.err"
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 env XMVB_PHASE_TIMING=1 "$source_dir/build/src/test_structure_tiles" --factor-generator 14 \
 > "$artifact/generator-phases.out" 2> "$artifact/generator-phases.err"
if [[ -n ${1:-} ]]; then
 bash "$source_dir/tests/mcscf/factor_gradient.sbatch" "$1"
fi
