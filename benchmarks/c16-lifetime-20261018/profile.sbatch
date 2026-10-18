#!/usr/bin/env bash
#SBATCH --job-name=C16-clock
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=00:10:00
#SBATCH --output=profile-%j.out
#SBATCH --error=profile-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0 cmake/3.31.9
set -u
native=$(realpath -e "$1")
artifact="$SLURM_SUBMIT_DIR/profile-$SLURM_JOB_ID"
mkdir -p "$artifact"
cd "$artifact"
cp "$native/native/C16_vb.xmi" "$native/native/cc-pvdz.gbs" .
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
export LD_LIBRARY_PATH=/home/guqqgroup/taoxia/xmvb-runs/porphine-current-20261017-c491230b/bundle/lib:${LD_LIBRARY_PATH:-}
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close
export OMP_MAX_ACTIVE_LEVELS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export XMVB_CIAH_TRACE=scalar
sha256sum "$native/bin/xmvb-cpp.exe" C16_vb.xmi cc-pvdz.gbs > hashes.txt
gprofng --version > version.txt 2>&1
set +e
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
  timeout -k 10s "${C16_PROFILE_SECONDS:-75}s" gprofng collect app -p on -o capture.er \
  stdbuf -oL -eL "$native/bin/xmvb-cpp.exe" C16_vb.xmi --optimizer-backend ciah \
  > output.out 2> run.err
printf '%s\n' "$?" > run.exit
set -e
gprofng display text -header -statistics capture.er > quality.txt 2> quality.err
gprofng display text -limit 40 -functions -calltree capture.er > hotspots.txt 2> report.err
