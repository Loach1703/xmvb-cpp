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
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|XMVB_PHASE_TIMING=|SLURM_CPUS_PER_TASK=)' > "$artifact/environment.txt"
cmake -S "$source_dir" -B "$source_dir/build" > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build" --target test_ci_factors test_structure_tiles test_orthogonal_direct_ci \
 test_coupled_structure test_product_frontier test_representation_overlap test_ci_transport \
 test_selected_ci_transport test_restricted_metric test_eigen_response test_keyframe -j8 > "$artifact/build.log" 2>&1
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
sha256sum "$source_dir/build/src/test_structure_tiles" > "$artifact/binary.sha256"
sha256sum "$source_dir/build/src/test_ci_factors" >> "$artifact/binary.sha256"
sha256sum "$source_dir/src/vbscf/structures/assembly/"{source_tiles.cpp,sources.cpp,sources.hpp,action.cpp} \
 "$source_dir/tests/vbscf/unit/test_structure_tiles.cpp" > "$artifact/sources.sha256"
sha256sum "$source_dir/src/vbscf/structures/assembly/"{source_response.cpp,source_direction.cpp,source_image.cpp,orthogonal_tiles.cpp,column_tiles.hpp,action.hpp} >> "$artifact/sources.sha256"
sha256sum "$source_dir/src/vbscf/derivatives/gradient/active_space/evaluator.cpp" \
 "$source_dir/src/vbscf/derivatives/hessian/exact/apply.cpp" \
 "$source_dir/src/vbscf/determinants/pairs/same_spin_cache.cpp" >> "$artifact/sources.sha256"
sha256sum "$source_dir/src/core/"{eigen_response.cpp,eigen_response.hpp} \
 "$source_dir/src/core/openmp.hpp" \
 "$source_dir/src/vbscf/structures/orthogonal_ci/exterior_transform.cpp" \
 "$source_dir/src/vbscf/derivatives/hessian/coupled/"{coupling.cpp,coupling.hpp} \
 "$source_dir/src/vbscf/derivatives/hessian/exact/"{operator.cpp,operator.hpp,state_internal.hpp} \
 "$source_dir/tests/vbscf/unit/test_keyframe.cpp" >> "$artifact/sources.sha256"
sha256sum "$source_dir/src/mcscf/ci/"{factors.cpp,factors.hpp,factors_internal.hpp,factor_generator.cpp,factor_plan.cpp,sigma.cpp,sigma_tiles.cpp,sigma.hpp,spin_channels.cpp} >> "$artifact/sources.sha256"
srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores env OMP_NUM_THREADS=2 \
 ctest --test-dir "$source_dir/build" \
 -R '^(ci_factors|orthogonal_direct_ci|coupled_structure|structure_tiles|product_frontier|representation_overlap|ci_transport|selected_ci_transport|restricted_metric|eigen_response|keyframe)$' \
 --output-on-failure > "$artifact/ctest.out" 2>&1
cp "$source_dir/build/Testing/Temporary/LastTest.log" "$artifact/test-details.log"
[[ ${1:-} == --tests-only ]] && exit 0
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 bash -c 'grep Cpus_allowed_list /proc/self/status' > "$artifact/affinity.txt"
if [[ ${1:-} == --accepted-state16 || ${1:-} == --accepted-hvp16 || ${1:-} == --profile-hvp16 ]]; then
 cmake --build "$source_dir/build" --target test_keyframe -j8 > "$artifact/accepted-build.log" 2>&1
 sha256sum "$source_dir/build/src/test_keyframe" >> "$artifact/binary.sha256"
 sha256sum "$source_dir/tests/vbscf/unit/test_keyframe.cpp" >> "$artifact/sources.sha256"
 if [[ $1 == --profile-hvp16 ]]; then
  # One bounded current-code capture. Profiled wall/RSS is not native timing.
  profiler=/opt/intel/oneapi/vtune/latest/bin64/vtune
  "$profiler" -version > "$artifact/profilers.txt" 2>&1
  set +e
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
   perf stat -e cycles:u,instructions:u -- true > "$artifact/perf-probe.out" 2> "$artifact/perf-probe.err"
  printf '%s\n' "$?" > "$artifact/perf-probe.exit"
  set -e
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
   env XMVB_PHASE_TIMING=1 "$profiler" -collect hotspots -knob sampling-mode=sw \
   -result-dir "$artifact/hvp.vtune" -- \
   "$source_dir/build/src/test_keyframe" --accepted-hvp16 \
   > "$artifact/profile.out" 2> "$artifact/profile.err"
  for report in summary hotspots callstacks; do
   srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
    "$profiler" -report "$report" -r "$artifact/hvp.vtune" \
    > "$artifact/profile-$report.txt" 2> "$artifact/profile-$report.err"
  done
  exit 0
 fi
 srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
  /usr/bin/time -v -o "$artifact/${1#--}.time" \
  "$source_dir/build/src/test_keyframe" "$1" \
  > "$artifact/${1#--}.out" 2> "$artifact/${1#--}.err"
 exit 0
fi
if [[ ${1:-} == --tile-adjoint12 ]]; then
 for mode in tile-adjoint12 tile-adjoint12-reference; do
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
   /usr/bin/time -v -o "$artifact/$mode.time" \
   "$source_dir/build/src/test_structure_tiles" --"$mode" \
   > "$artifact/$mode.out" 2> "$artifact/$mode.err"
 done
 exit 0
fi
if [[ ${1:-} == --production || ${1:-} == --selected-scf ]]; then
 if [[ ${1:-} == --production ]]; then
  bash "$source_dir/tests/mcscf/factor_gradient.sbatch" "$2"
 else
  cmake --build "$source_dir/build" --target xmvb -j8 > "$artifact/application-build.log" 2>&1
 fi
 sha256sum "$source_dir/build/src/xmvb-cpp.exe" >> "$artifact/binary.sha256"
 python3 "$source_dir/tests/mcscf/selected_scf.py" prepare "$2" "$artifact" --system "${3:-240}"
 for count in 1 3; do
  for solver in davidson dense; do
   cd "$artifact/selected$count-$solver"
   sha256sum input.xmi cc-pvdz.gbs > input.sha256
   srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
    /usr/bin/time -v -o resource.time "$source_dir/build/src/xmvb-cpp.exe" input.xmi \
    --optimizer-backend neo --eigensolver "$solver" --orbital-preconditioner hessian-diagonal \
    --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --neo-trace neo.tsv \
    > output.out 2> run.err
  done
 done
 python3 "$source_dir/tests/mcscf/selected_scf.py" summarize "$artifact" "$artifact/selected_scf.csv"
 exit 0
fi
if [[ ${1:-} == --source-states ]]; then
 for n in 12 14; do
  for mode in source-state source-state-reference; do
   srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
    /usr/bin/time -v -o "$artifact/$n-$mode.time" \
    "$source_dir/build/src/test_structure_tiles" --"$mode" "$n" \
    > "$artifact/$n-$mode.out" 2> "$artifact/$n-$mode.err"
  done
  count=1
  [[ $n == 14 ]] && count=2
  for mode in structure-response$n structure-response$n-reference; do
   srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
    /usr/bin/time -v -o "$artifact/$count-$mode.time" \
    "$source_dir/build/src/test_structure_tiles" --"$mode" "$count" \
    > "$artifact/$count-$mode.out" 2> "$artifact/$count-$mode.err"
  done
 done
 srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
  /usr/bin/time -v -o "$artifact/16-source-state.time" \
  "$source_dir/build/src/test_structure_tiles" --source-accepted16 \
  > "$artifact/16-source-state.out" 2> "$artifact/16-source-state.err"
 if [[ -n ${2:-} ]]; then
  bash "$source_dir/tests/mcscf/factor_gradient.sbatch" "$2"
 fi
 exit 0
fi
if [[ ${1:-} == --response14 ]]; then
 for mode in structure-response14 structure-response14-reference; do
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
   /usr/bin/time -v -o "$artifact/2-$mode.time" \
   "$source_dir/build/src/test_structure_tiles" --"$mode" 2 \
   > "$artifact/2-$mode.out" 2> "$artifact/2-$mode.err"
 done
 exit 0
fi
if [[ ${1:-} == --responses ]]; then
 for count in 1 3 20 132; do
  for mode in structure-response12 structure-response12-reference; do
   srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
    /usr/bin/time -v -o "$artifact/$count-$mode.time" \
    "$source_dir/build/src/test_structure_tiles" --"$mode" "$count" \
    > "$artifact/$count-$mode.out" 2> "$artifact/$count-$mode.err"
  done
 done
 if [[ -n ${2:-} ]]; then
  bash "$source_dir/tests/mcscf/factor_gradient.sbatch" "$2"
 fi
 exit 0
fi
if [[ ${1:-} == --product-response14 ]]; then
 for suffix in '' '-reference'; do
  srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
   /usr/bin/time -v -o "$artifact/132-structure-response14$suffix.time" \
   "$source_dir/build/src/test_structure_tiles" --structure-response14-once"$suffix" 132 \
   > "$artifact/132-structure-response14$suffix.out" 2> "$artifact/132-structure-response14$suffix.err"
 done
 exit 0
fi
if [[ ${1:-} == --generator16 ]]; then
 srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
  /usr/bin/time -v -o "$artifact/16-factor-generator.time" \
  "$source_dir/build/src/test_structure_tiles" --factor-generator 16 \
  > "$artifact/16-factor-generator.out" 2> "$artifact/16-factor-generator.err"
 exit 0
fi
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
 /usr/bin/time -v -o "$artifact/16-factor-generator.time" \
 "$source_dir/build/src/test_structure_tiles" --factor-generator 16 \
 > "$artifact/16-factor-generator.out" 2> "$artifact/16-factor-generator.err"
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 env XMVB_PHASE_TIMING=1 /usr/bin/time -v -o "$artifact/16-direction-check.time" \
 "$source_dir/build/src/test_structure_tiles" --factor-generator-check16 \
 > "$artifact/16-direction-check.out" 2> "$artifact/16-direction-check.err"
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
