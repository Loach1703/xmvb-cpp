#!/usr/bin/env bash
#SBATCH --job-name=casscf-coarse
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --array=0-2
#SBATCH --output=run-%A_%a.out
#SBATCH --error=run-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0
set -u
source_dir=$(realpath -e "$1")
old_runs="$source_dir/benchmarks/casscf-20261017"
build="$source_dir/build-final/src"
bundle=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle
export LD_LIBRARY_PATH="$bundle/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
# Leave the normal run uninstrumented; explicit scalar tracing is available
# for a separate diagnostic allocation, never for the native timing claim.
case $SLURM_ARRAY_TASK_ID in
  0) case_name=native-F2; input="$old_runs/native-F2-291885/input.xmi" ;;
  1) case_name=native-FeCl2; input="$old_runs/native-FeCl2-291886/input.xmi" ;;
  2) case_name=fixture-FeCl2; input="$old_runs/casscf-FeCl2-291883/fixture/ao-mo.bin" ;;
  3) case_name=vb-F2; input="$source_dir/testdata/vbscf/F2.xmi" ;;
  4) case_name=vb-FeCl2; input="$source_dir/testdata/vbscf/FeCl2.xmi" ;;
  *) exit 64 ;;
esac
artifact="$SLURM_SUBMIT_DIR/${case_name}-${SLURM_JOB_ID}"
mkdir "$artifact"
cd "$artifact"
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
cp "$source_dir/source.sha256" .
if [[ $case_name == fixture-* ]]; then
  command=("$build/casscf_fixture" "$input" "$artifact/trajectory.csv")
else
  cp "$input" input.xmi
  command=("$build/xmvb-cpp.exe" input.xmi)
  if [[ $case_name == vb-* ]]; then
    command+=(--optimizer-backend ciah --eigensolver davidson
      --orbital-preconditioner hessian-diagonal)
  fi
fi
sha256sum "${command[0]}" "$input" > hashes.txt
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  bash -c 'grep Cpus_allowed_list /proc/self/status' > affinity.txt
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time "${command[@]}" > output.out 2> run.err
if [[ $case_name == fixture-* ]]; then
  grep -q '^converged=1$' output.out
elif [[ $case_name == vb-* ]]; then
  grep -q '^Termination reason *: ciah_projected_gradient_tolerance' output.out
else
  grep -q '^Converged: yes$' output.out
  [[ -s input.orb ]]
fi
