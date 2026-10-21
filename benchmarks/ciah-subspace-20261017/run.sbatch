#!/usr/bin/env bash
#SBATCH --job-name=ciah-subspace-singlet
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=192G
#SBATCH --time=15-00:00:00
#SBATCH --output=singlet-%j.out
#SBATCH --error=singlet-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 && $SLURM_CPUS_PER_TASK == 32 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0
set -u
source_dir=$(realpath -e "$1")
baseline=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle
run="$SLURM_SUBMIT_DIR/singlet-$SLURM_JOB_ID"
mkdir -p "$run"
cp "$baseline/cases/singlet/input.xmi" "$baseline/cc-pvdz.gbs" "$run/"
cp "$source_dir/source.sha256" "$source_dir/parent-revision.txt" "$run/"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="$baseline/lib:${LD_LIBRARY_PATH:-}"
# Scalar events use existing images only, unlike the expensive audit mode.
export XMVB_CIAH_TRACE=scalar
unset XMVB_NEO_TRACE XMVB_NEO_KEYFRAME_AUDIT XMVB_NEO_ORBITAL_STEP_AUDIT
unset XMVB_CIAH_CANDIDATE_AUDIT XMVB_NEO_CANDIDATE_AUDIT XMVB_PHASE_TIMING
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_|XMVB_)' > environment.txt
sha256sum "$source_dir/build/src/xmvb-cpp.exe" input.xmi cc-pvdz.gbs "$baseline/lib/"* > hashes.txt
arguments=(--optimizer-backend ciah --eigensolver davidson
  --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3
  --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 --ciah-trace ciah.tsv)
printf '%s\n' "${arguments[@]}" > arguments.txt
export CIAH_SUBSPACE_EXECUTABLE="$source_dir/build/src/xmvb-cpp.exe"
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -c '
  grep Cpus_allowed_list /proc/self/status > affinity.txt
  exec /usr/bin/time -v -o resource.time "$CIAH_SUBSPACE_EXECUTABLE" input.xmi "$@"
' bash "${arguments[@]}" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
[[ $status == 0 ]] || exit "$status"
grep -q '^Termination reason *: ciah_projected_gradient_tolerance' output.out
