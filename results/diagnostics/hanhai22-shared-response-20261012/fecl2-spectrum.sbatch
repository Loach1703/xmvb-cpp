#!/usr/bin/env bash
#SBATCH --job-name=fecl2-spectrum
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --output=fecl2-spectrum-%j.out
#SBATCH --error=fecl2-spectrum-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
set -u
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
run="$root/fecl2-spectrum-$SLURM_JOB_ID"
mkdir -p "$run"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > "$run/node.txt"
lscpu > "$run/cpu.txt"
cmake --build "$source_dir/build" --target xmvb test_ciah test_coupled_structure test_orthogonal_direct_ci test_sparse_orbital_quotient test_hessian_diagonal -j8 > "$run/build.log" 2>&1
for test_name in test_ciah test_coupled_structure test_orthogonal_direct_ci test_sparse_orbital_quotient test_hessian_diagonal; do
 srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores "$source_dir/build/src/$test_name" > "$run/$test_name.out" 2> "$run/$test_name.err"
done
cp "$source_dir/build/src/xmvb-cpp.exe" "$run/"
cp "$root/fecl2-ritz-cost-1248947/native/input.xmi" "$run/"
cp "$root/fecl2-ritz-cost-1248947/native/def2-tzvp.gbs" "$run/"
export XMVB_NEO_TRACE=1 XMVB_NEO_CANDIDATE_AUDIT=1
cd "$run"
sha256sum xmvb-cpp.exe input.xmi def2-tzvp.gbs > hashes.txt
status=0
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
 ./xmvb-cpp.exe input.xmi --optimizer-backend neo --eigensolver davidson \
 --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
 --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
[[ "$status" == 0 || "$status" == 2 ]]
grep -q 'ciah_orbital_preconditioner_spectrum' run.err
