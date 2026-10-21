#!/usr/bin/env bash
#SBATCH --job-name=porphine-current
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=192G
#SBATCH --time=15-00:00:00
#SBATCH --output=porphine-current-%A_%a.out
#SBATCH --error=porphine-current-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 && $SLURM_CPUS_PER_TASK == 32 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
set -u
bundle=$(realpath -e "$1")
(cd "$bundle" && sha256sum --check hashes.sha256) > /dev/null
states=(singlet triplet quintet)
state=${states[$SLURM_ARRAY_TASK_ID]}
run="$SLURM_SUBMIT_DIR/porphine-current-$SLURM_ARRAY_JOB_ID/$state"
mkdir -p "$run"
cp "$bundle/cases/$state/input.xmi" "$bundle/cc-pvdz.gbs" "$run/"
cp "$bundle/"{revision.txt,sources.sha256,ctest.out} "$run/"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE XMVB_NEO_TRACE XMVB_NEO_KEYFRAME_AUDIT
unset XMVB_NEO_ORBITAL_STEP_AUDIT XMVB_NEO_ORBITAL_PATH XMVB_NEO_CANDIDATE_AUDIT
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
printf '%s\n' "$bundle" > bundle.txt
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > environment.txt
sha256sum "$bundle/xmvb-cpp.exe" input.xmi cc-pvdz.gbs > hashes.txt
arguments=(--optimizer-backend ciah --eigensolver davidson
  --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3
  --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 --ciah-trace ciah.tsv)
printf '%s\n' "${arguments[@]}" > arguments.txt
export PORPHINE_CURRENT_BINARY="$bundle/xmvb-cpp.exe"
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -c '
  grep Cpus_allowed_list /proc/self/status > affinity.txt
  exec /usr/bin/time -v -o resource.time "$PORPHINE_CURRENT_BINARY" input.xmi "$@"
' bash "${arguments[@]}" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
[[ $status == 0 ]] || exit "$status"
grep -q '^Termination reason *: ciah_projected_gradient_tolerance' output.out
