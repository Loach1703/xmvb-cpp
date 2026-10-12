#!/usr/bin/env bash
#SBATCH --job-name=fecl2-ritz-cost
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=fecl2-ritz-cost-%j.out
#SBATCH --error=fecl2-ritz-cost-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
artifact="$root/fecl2-ritz-cost-$SLURM_JOB_ID"
mkdir -p "$artifact/native" "$artifact/profile"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$artifact/node.txt"
lscpu > "$artifact/cpu.txt"
module list > "$artifact/modules.txt" 2>&1
cmake --build "$source_dir/build" --target xmvb test_ciah test_coupled_structure test_orthogonal_direct_ci -j8 > "$artifact/build.log" 2>&1
for test_name in test_ciah test_coupled_structure test_orthogonal_direct_ci; do
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores "$source_dir/build/src/$test_name" > "$artifact/$test_name.out" 2> "$artifact/$test_name.err"
done
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/"
grep -E 'CMAKE_BUILD_TYPE:|CMAKE_CXX_FLAGS(_RELEASE)?:' "$source_dir/build/CMakeCache.txt" > "$artifact/flags.txt"
for variant in native profile; do
 cp "$root/fecl2-metric-convergence-1248707/input.xmi" "$artifact/$variant/"
 cp "$root/fecl2-metric-convergence-1248707/def2-tzvp.gbs" "$artifact/$variant/"
done
cd "$artifact/native"
sha256sum ../xmvb-cpp.exe input.xmi def2-tzvp.gbs > hashes.txt
unset XMVB_NEO_CANDIDATE_AUDIT
export XMVB_NEO_TRACE=scalar
status=0
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
 ../xmvb-cpp.exe input.xmi --optimizer-backend neo --eigensolver davidson \
 --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
 --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
[[ "$status" == 0 || "$status" == 2 ]]
[[ ${XMVB_AUDIT_SKIP_PROFILE:-0} == 1 ]] && exit 0
cd "$artifact/profile"
# Bound the independent software sampling capture to eight original-guess macros.
sed -i -E 's/itmax=[0-9]+/itmax=8/g' input.xmi
unset XMVB_NEO_TRACE
status=0
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
 /opt/intel/oneapi/vtune/latest/bin64/vtune -collect hotspots -knob sampling-mode=sw \
 -result-dir capture -- ../xmvb-cpp.exe input.xmi --optimizer-backend neo --eigensolver davidson \
 --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
 --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
if [[ -d capture ]]; then
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
 /opt/intel/oneapi/vtune/latest/bin64/vtune -report hotspots -r capture > hotspots.txt 2> report.err
fi
