#!/usr/bin/env bash
#SBATCH --job-name=porphine-keyframe-audit
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=01:00:00
#SBATCH --output=porphine-keyframe-audit-%j.out
#SBATCH --error=porphine-keyframe-audit-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
run="$root/porphine-keyframe-audit-$SLURM_JOB_ID"
mkdir -p "$run"
cd "$run"
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
hostname > node.txt
sha256sum "$source_dir/src/vbscf/optimization/backends/ciah.cpp" \
  "$source_dir/src/vbscf/diagnostics/orbitals/ciah_keyframe.cpp" \
  "$source_dir/src/vbscf/derivatives/hessian/coupled/structure.cpp" > source-hashes.txt
cmake -S "$source_dir" -B "$source_dir/build" > configure.log 2>&1
cmake --build "$source_dir/build" --target xmvb test_ciah test_coupled_structure \
  -j16 > build.log 2>&1
for test in test_ciah test_coupled_structure; do
  srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
    "$source_dir/build/src/$test" > "$test.out" 2> "$test.err"
done
cp "$source_dir/build/src/xmvb-cpp.exe" diagnostic.exe
cp "$root/inputs/singlet/input.xmi" "$root/inputs/singlet/cc-pvdz.gbs" .
sha256sum diagnostic.exe input.xmi cc-pvdz.gbs > hashes.txt
export OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export XMVB_NEO_TRACE=scalar XMVB_NEO_KEYFRAME_AUDIT=1
unset XMVB_NEO_CANDIDATE_AUDIT
status=0
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  /usr/bin/time -v -o resource.time ./diagnostic.exe input.xmi \
  --optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal \
  --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 \
  --max-iterations 1 --neo-trace neo.tsv > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
if [[ "$status" == 2 ]]; then exit 0; fi
exit "$status"
