#!/usr/bin/env bash
#SBATCH --job-name=C14-metric-fix
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=128G
#SBATCH --time=15-00:00:00
#SBATCH --output=C14-metric-fix-%j.out
#SBATCH --error=C14-metric-fix-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
set -u
source_dir=/home/guqqgroup/xiatao/work/xmvb-source-tiles-20261015
artifact="$SLURM_SUBMIT_DIR/fix-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
cmake -S "$source_dir" -B "$source_dir/build" > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build" --target xmvb test_selected_metric \
 test_selected_ci_transport test_coupled_structure test_representation_overlap \
 test_ci_transport test_keyframe test_ciah -j16 > "$artifact/build.log" 2>&1
cp "$source_dir/build/CMakeCache.txt" "$artifact/"
srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores env OMP_NUM_THREADS=2 \
 ctest --test-dir "$source_dir/build" --output-on-failure \
 -R '^(selected_metric|selected_ci_transport|coupled_structure|representation_overlap|ci_transport|keyframe|ciah)$' \
 > "$artifact/ctest.out" 2>&1
cp "$source_dir/build/Testing/Temporary/LastTest.log" "$artifact/test-details.log"
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/xmvb-cpp.exe"
cp "$SLURM_SUBMIT_DIR/C14.xmi" "$artifact/C14.xmi"
ln -s "$source_dir/basis" "$artifact/basis"
cd "$artifact"
sha256sum xmvb-cpp.exe C14.xmi basis/cc-pvdz.gbs > hashes.sha256
env | sort | sed -n '/^OMP_/p; /^OPENBLAS_/p; /^MKL_/p; /^BLIS_/p; /^SLURM_/p' > environment.txt
find "$source_dir/src" "$source_dir/tests/vbscf" -type f \
 \( -name '*.cpp' -o -name '*.hpp' -o -name '*.cmake' \) -print0 \
 | sort -z | xargs -0 sha256sum > sources.sha256
set +e
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -c '
  grep Cpus_allowed_list /proc/self/status > affinity.txt
  exec /usr/bin/time -v -o resource.time \
    ./xmvb-cpp.exe C14.xmi --ciah-trace ciah.tsv
' > C14.out 2> run.err
status=$?
printf '%s\n' "$status" > run.exit
set -e
[[ "$status" == 0 ]]
grep -q '^Termination reason *: ciah_projected_gradient_tolerance' C14.out
