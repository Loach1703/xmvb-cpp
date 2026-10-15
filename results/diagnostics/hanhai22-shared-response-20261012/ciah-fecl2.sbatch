#!/usr/bin/env bash
#SBATCH --job-name=fecl2-ciah
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=ciah-%j.out
#SBATCH --error=ciah-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
artifact="$root/ciah-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
cmake -S "$source_dir" -B "$source_dir/build" > "$artifact/configure.log" 2>&1
cmake --build "$source_dir/build" --target xmvb test_ciah test_neo_response_solver test_coupled_structure -j8 > "$artifact/build.log" 2>&1
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/"
cp "$root/fecl2-first-step-1247319/input.xmi" "$root/fecl2-first-step-1247319/def2-tzvp.gbs" "$artifact/"
cd "$artifact"
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
cp "$source_dir/build/src/CMakeFiles/xmvb.dir/flags.make" build-flags.txt 2>/dev/null || true
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > execution-env.txt
for test_name in test_ciah test_neo_response_solver test_coupled_structure; do
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
   "$source_dir/build/src/$test_name" > "$test_name.out" 2> "$test_name.err"
done
sha256sum xmvb-cpp.exe input.xmi def2-tzvp.gbs > hashes.txt
set +e
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
 env XMVB_NEO_TRACE=scalar /usr/bin/time -v -o resource.time \
 ./xmvb-cpp.exe input.xmi --neo-trace neo.tsv > output.out 2> run.err
printf '%s\n' "$?" > run.exit
