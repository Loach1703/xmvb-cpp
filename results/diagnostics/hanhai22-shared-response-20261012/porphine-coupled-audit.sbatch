#!/usr/bin/env bash
#SBATCH --job-name=porphine-coupled-audit
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=01:00:00
#SBATCH --output=porphine-coupled-audit-%j.out
#SBATCH --error=porphine-coupled-audit-%j.err
set -eo pipefail
. /etc/profile
module purge
module load gcc/14.1.0 cmake/3.30.5
root="$SLURM_SUBMIT_DIR"
source_dir=/home/guqqgroup/xiatao/work/xmvb-shared-response-1eb1b50a
artifact="$root/porphine-coupled-audit-$SLURM_JOB_ID"
mkdir -p "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
cmake --build "$source_dir/build" --target audit_orbital_gradient xmvb -j8 > "$artifact/build.log" 2>&1
cd "$artifact"
cp "$root/porphine-space-gradient-1247942/endpoint.xmi" .
cp "$root/inputs/triplet/cc-pvdz.gbs" .
sha256sum "$source_dir/build/src/audit_orbital_gradient" "$source_dir/build/src/xmvb-cpp.exe" endpoint.xmi > hashes.txt
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o coupled.time \
 "$source_dir/build/src/audit_orbital_gradient" endpoint.xmi davidson 1e-7 1e-3 gradient coupled-keyframe \
 > coupled.out 2> coupled.err
# This is a bounded restart from the exported accepted point. It does not
# recreate the original failed macrostep's unsaved seed or bit-exact orbitals.
python3 -c 'import pathlib; p=pathlib.Path("endpoint.xmi"); pathlib.Path("restart.xmi").write_text(p.read_text().replace("itmax=1000", "itmax=2"))'
export XMVB_NEO_TRACE=scalar
status=0
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores /usr/bin/time -v -o restart.time \
 "$source_dir/build/src/xmvb-cpp.exe" restart.xmi --optimizer-backend neo --eigensolver davidson \
 --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
 --exact-exchange-memory-mib 256 --neo-trace restart.tsv > restart.out 2> restart.err || status=$?
printf '%s\n' "$status" > restart.exit
