#!/usr/bin/env bash
#SBATCH --job-name=simd-ri-pair
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=15-00:00:00
#SBATCH --output=simd-ri-%j.out
#SBATCH --error=simd-ri-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
baseline="$1"
updated="$2"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
root="$SLURM_SUBMIT_DIR/simd-ri-$SLURM_JOB_ID"
mkdir -p "$root"
hostname > "$root/node.txt"
for variant in baseline updated; do
    source_dir="$baseline"
    [[ "$variant" != updated ]] || source_dir="$updated"
    mkdir -p "$root/$variant"
    cd "$root/$variant"
    cp "$baseline/testdata/vbscf/F2.xmi" input.xmi
    # Mechanical change only in this disposable run input, not the source deck.
    sed -i 's/INT=LIBCINT/INT=RI/' input.xmi
    exe="$source_dir/build/src/xmvb-cpp.exe"
    sha256sum "$exe" input.xmi > hashes.txt
    srun --exact --cpu-bind=cores /usr/bin/time -v -o resource.time "$exe" input.xmi \
        --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        > output.out 2> run.err
done
