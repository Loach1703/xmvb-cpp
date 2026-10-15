#!/usr/bin/env bash
#SBATCH --job-name=simd-shell-heap
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=15-00:00:00
#SBATCH --output=simd-heap-%j.out
#SBATCH --error=simd-heap-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
baseline="$1"
updated="$2"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export PATH="/home/guqqgroup/taoxia/software/bin:$dependencies/bin:$PATH"
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=8 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
root="$SLURM_SUBMIT_DIR/simd-heap-$SLURM_JOB_ID"
mkdir -p "$root"
hostname > "$root/node.txt"
for variant in baseline updated; do
    source_dir="$baseline"
    [[ "$variant" != updated ]] || source_dir="$updated"
    mkdir -p "$root/$variant"
    cd "$root/$variant"
    cp "$baseline/testdata/vbscf/241_vb.xmi" input.xmi
    exe="$source_dir/build/src/xmvb-cpp.exe"
    sha256sum "$exe" input.xmi > hashes.txt
    srun --exact --cpu-bind=cores heaptrack -o memory "$exe" input.xmi \
        --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        > output.out 2> run.err
    shopt -s nullglob
    captures=(memory*.gz memory*.zst)
    [[ ${#captures[@]} == 1 ]] || exit 65
    srun --exact --cpu-bind=cores heaptrack_print "${captures[0]}" \
        > memory-report.txt 2> analysis.err
done
