#!/usr/bin/env bash
#SBATCH --job-name=240-metric
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=240-metric-%j.out
#SBATCH --error=240-metric-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && -d ${1:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
source_bundle=$(realpath "$1")
shift
export AUDIT_ROOT="$SLURM_SUBMIT_DIR/240-metric-$SLURM_JOB_ID"
mkdir -p "$AUDIT_ROOT/bundle"
cmake -S "$SLURM_SUBMIT_DIR" -B "$SLURM_SUBMIT_DIR/build" > "$AUDIT_ROOT/configure.log" 2>&1
cmake --build "$SLURM_SUBMIT_DIR/build" --target xmvb test_product_frontier \
 test_structure_tiles test_representation_overlap test_ci_transport test_selected_ci_transport \
 "$@" -j8 > "$AUDIT_ROOT/build.log" 2>&1
cp "$SLURM_SUBMIT_DIR/build/src/xmvb-cpp.exe" "$source_bundle/240_vb.xmi" \
 "$source_bundle/cc-pvdz.gbs" "$AUDIT_ROOT/bundle/"
cp "$SLURM_SUBMIT_DIR/build/CMakeCache.txt" "$AUDIT_ROOT/"
sha256sum "$AUDIT_ROOT/bundle/xmvb-cpp.exe" > "$AUDIT_ROOT/binary.sha256"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
hostname > "$AUDIT_ROOT/node.txt"
lscpu > "$AUDIT_ROOT/cpu.txt"
srun --exact --ntasks=1 --cpus-per-task=2 --cpu-bind=cores env OMP_NUM_THREADS=2 \
 ctest --test-dir "$SLURM_SUBMIT_DIR/build" \
 -R '^(product_frontier|structure_tiles|representation_overlap|ci_transport|selected_ci_transport|selected_metric)$' \
 --output-on-failure > "$AUDIT_ROOT/ctest.out" 2>&1
srun --exact --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores bash -se <<'RUN'
grep Cpus_allowed_list /proc/self/status > "$AUDIT_ROOT/affinity.txt"
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_|XMVB_)' > "$AUDIT_ROOT/environment.txt"
labels=(cov ion01 ion02 ion03 full)
keywords=('COV' 'ION(0:1)' 'ION(0:2)' 'ION(0:3)' 'FULL')
for mode in native phases; do
 for index in 0 1 2 3 4; do
  [[ $mode == phases && $index -lt 2 ]] && continue
  mkdir -p "$AUDIT_ROOT/$mode/${labels[index]}"
  cd "$AUDIT_ROOT/$mode/${labels[index]}"
  sed -E -e "s/^str=.*/str=${keywords[index]}/I" -e 's/^iscf=.*/iscf=9/I' \
   "$AUDIT_ROOT/bundle/240_vb.xmi" > input.xmi
  cp "$AUDIT_ROOT/bundle/cc-pvdz.gbs" .
  unset XMVB_PHASE_TIMING
  [[ $mode == phases ]] && export XMVB_PHASE_TIMING=1
  sha256sum input.xmi cc-pvdz.gbs > input.sha256
  /usr/bin/time -v -o resource.time "$AUDIT_ROOT/bundle/xmvb-cpp.exe" input.xmi \
   --optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal \
   --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --neo-trace neo.tsv > output.out 2> run.err
  printf '0\n' > run.exit
 done
done
RUN
