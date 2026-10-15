#!/usr/bin/env bash
#SBATCH --job-name=fecl2-independent-seed
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=fecl2-seed-%j.out
#SBATCH --error=fecl2-seed-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
artifact="$root/fecl2-seed-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
cmake --build "$source_dir/build" --target xmvb test_neo_response_solver -j8 > "$artifact/build.log" 2>&1
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/"
cp "$root/fecl2-first-step-1247319/input.xmi" "$root/fecl2-first-step-1247319/def2-tzvp.gbs" "$artifact/"
cd "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
set +e
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  "$source_dir/build/src/test_neo_response_solver" > unit.out 2> unit.err
printf '%s\n' "$?" > unit.exit
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  env XMVB_NEO_TRACE=scalar /usr/bin/time -v -o resource.time \
  ./xmvb-cpp.exe input.xmi --max-iterations 30 --neo-trace neo.tsv \
  > output.out 2> run.err
printf '%s\n' "$?" > run.exit
