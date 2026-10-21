#!/usr/bin/env bash
#SBATCH --job-name=RI-retained-build
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --mem=96G
#SBATCH --time=00:10:00
#SBATCH --output=retained-build-%j.out
#SBATCH --error=retained-build-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} && $# == 1 ]] || exit 64
. /etc/profile
module purge
module load gcc/14.3.0 cmake/3.31.9
set -u
source_dir=$(realpath -e "$1")
artifact="$SLURM_SUBMIT_DIR/retained-build-$SLURM_JOB_ID"
mkdir -p "$artifact"
cd "$source_dir"
cp "$SLURM_SUBMIT_DIR/retained-source.sha256" "$artifact/source.sha256"
sha256sum --strict --check "$artifact/source.sha256" > "$artifact/source-check.txt"
cmake --build build-final --target xmvb -j"$SLURM_CPUS_PER_TASK" > "$artifact/build.log" 2>&1
sha256sum build-final/src/xmvb-cpp.exe > "$artifact/binary.sha256"
