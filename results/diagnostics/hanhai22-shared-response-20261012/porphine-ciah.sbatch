#!/usr/bin/env bash
#SBATCH --job-name=porphine-ciah-diag
#SBATCH --partition=CPU-192C768GB
#SBATCH --qos=qos_cpu_192c768gb
#SBATCH --account=xiatao
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=15-00:00:00
#SBATCH --output=porphine-ciah-%A_%a.out
#SBATCH --error=porphine-ciah-%A_%a.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load gcc/14.1.0
root="$SLURM_SUBMIT_DIR"
states=(singlet triplet)
state="${states[$SLURM_ARRAY_TASK_ID]}"
frozen="$root/ciah-positive-build-1247629"
binary="$frozen/xmvb-cpp.exe"
expected=10430688c99848326f6c21967d4bf7e43c4a2ba1db1628b6e7a63b1e064e98a1
actual=$(sha256sum "$binary")
[[ "${actual%% *}" == "$expected" ]] || exit 65
artifact="$root/porphine-ciah-$SLURM_ARRAY_JOB_ID/$state"
mkdir -p "$artifact"
cp "$root/inputs/$state/input.xmi" "$root/inputs/$state/cc-pvdz.gbs" "$artifact/"
cp "$frozen/build-flags.txt" "$frozen/source.txt" "$artifact/"
cd "$artifact"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/guqqgroup/xiatao/software/xmvb-deps-20261012/lib:${LD_LIBRARY_PATH:-}"
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
env | sort | grep -E '^(OMP_|OPENBLAS_|MKL_|BLIS_|SLURM_)' > execution-env.txt
sha256sum "$binary" input.xmi cc-pvdz.gbs > hashes.txt
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
 "$root/affinity-1246753/affinity_probe" > affinity.out 2> affinity.err
printf '%s\n' "$binary input.xmi --optimizer-backend neo --eigensolver davidson --orbital-preconditioner hessian-diagonal --gradient-tolerance 1e-3 --energy-tolerance 1e-7 --exact-exchange-memory-mib 256 --neo-trace neo.tsv" > command.txt
date -Is > started.txt
set +e
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
 /usr/bin/time -v -o resource.time "$binary" input.xmi \
 --optimizer-backend neo --eigensolver davidson \
 --orbital-preconditioner hessian-diagonal \
 --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
 --exact-exchange-memory-mib 256 --neo-trace neo.tsv > output.out 2> run.err
status=$?
printf '%s\n' "$status" > run.exit
date -Is > finished.txt
exit "$status"
