#!/usr/bin/env bash
#SBATCH --job-name=ciah-macro-audit
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=01:00:00
#SBATCH --output=ciah-macro-audit-%j.out
#SBATCH --error=ciah-macro-audit-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
artifact="$root/ciah-macro-audit-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
cmake --build "$source_dir/build" --target xmvb test_ciah -j8 > "$artifact/build.log" 2>&1
"$source_dir/build/src/test_ciah" > "$artifact/unit.out" 2> "$artifact/unit.err"
cp "$source_dir/build/src/xmvb-cpp.exe" "$artifact/"
sha256sum "$artifact/xmvb-cpp.exe" > "$artifact/hash.txt"
cd "$artifact"
cp "$root/inputs/triplet/cc-pvdz.gbs" .
python3 -c 'import pathlib; p=pathlib.Path("'"$root"'/porphine-space-gradient-1247942/endpoint.xmi"); pathlib.Path("restart.xmi").write_text(p.read_text().replace("itmax=1000", "itmax=2"))'
export XMVB_NEO_TRACE=1 XMVB_NEO_CANDIDATE_AUDIT=1
status=0
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o resource.time \
 "$artifact/xmvb-cpp.exe" restart.xmi --optimizer-backend neo --eigensolver davidson \
 --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
 --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err || status=$?
printf '%s\n' "$status" > run.exit
# Native exit 2 is the intentional diagnostic iteration limit, not convergence.
if [[ "$status" != 0 && "$status" != 2 ]]; then exit "$status"; fi
