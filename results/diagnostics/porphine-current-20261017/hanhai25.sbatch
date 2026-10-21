#!/usr/bin/env bash
#SBATCH --job-name=porphine-current
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=192G
#SBATCH --time=15-00:00:00
#SBATCH --output=porphine-h25-%A_%a.out
#SBATCH --error=porphine-h25-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# -ge 1 && $SLURM_CPUS_PER_TASK == 32 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0
set -u
bundle=$(realpath -e "$1")
(cd "$bundle" && sha256sum --check hashes-h25.sha256) > /dev/null
if [[ ${2:-} == --qualify ]]; then
  state=qualification
  run="$SLURM_SUBMIT_DIR/qualification-$SLURM_JOB_ID"
else
  states=(singlet triplet quintet)
  state=${states[$SLURM_ARRAY_TASK_ID]}
  run="$SLURM_SUBMIT_DIR/run-$SLURM_ARRAY_JOB_ID/$state"
fi
mkdir -p "$run"
cp "$bundle/cc-pvdz.gbs" "$run/"
cp "$bundle/"{revision.txt,sources.sha256,ctest.out,hashes-h25.sha256} "$run/"
cd "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="$bundle/lib:${LD_LIBRARY_PATH:-}"
unset XMVB_PHASE_TIMING XMVB_CIAH_TRACE XMVB_NEO_TRACE XMVB_NEO_KEYFRAME_AUDIT
unset XMVB_NEO_ORBITAL_STEP_AUDIT XMVB_NEO_ORBITAL_PATH XMVB_NEO_CANDIDATE_AUDIT
export PORPHINE_CURRENT_BUNDLE="$bundle"
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
printf '%s\n' "$bundle" > bundle.txt
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > environment.txt
arguments=(--optimizer-backend ciah --eigensolver davidson
  --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3
  --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 --ciah-trace ciah.tsv)
printf '%s\n' "${arguments[@]}" > arguments.txt
if [[ $state == qualification ]]; then
  srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -se <<'CHECK'
grep Cpus_allowed_list /proc/self/status > affinity.txt
ldd "$PORPHINE_CURRENT_BUNDLE/xmvb-cpp.exe" > dependencies.txt
! grep -q 'not found' dependencies.txt
env OMP_NUM_THREADS=2 "$PORPHINE_CURRENT_BUNDLE/test_structure_tiles" > tiles.out 2> tiles.err
cp "$PORPHINE_CURRENT_BUNDLE/smoke/input.xmi" .
/usr/bin/time -v -o resource.time "$PORPHINE_CURRENT_BUNDLE/xmvb-cpp.exe" input.xmi \
  --optimizer-backend ciah --eigensolver davidson --orbital-preconditioner hessian-diagonal \
  --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 \
  --ciah-trace ciah.tsv > output.out 2> run.err
grep -q '^Termination reason *: ciah_projected_gradient_tolerance' output.out
CHECK
  exit 0
fi
cp "$bundle/cases/$state/input.xmi" .
sha256sum "$bundle/xmvb-cpp.exe" input.xmi cc-pvdz.gbs "$bundle/lib/"* > hashes.txt
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores bash -c '
  grep Cpus_allowed_list /proc/self/status > affinity.txt
  exec /usr/bin/time -v -o resource.time "$PORPHINE_CURRENT_BUNDLE/xmvb-cpp.exe" input.xmi "$@"
' bash "${arguments[@]}" > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
[[ $status == 0 ]] || exit "$status"
grep -q '^Termination reason *: ciah_projected_gradient_tolerance' output.out
