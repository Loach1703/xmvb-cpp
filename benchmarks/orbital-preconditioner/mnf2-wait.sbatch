#!/usr/bin/env bash
#SBATCH --job-name=mnf2-wait
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15-00:00:00
#SBATCH --output=mnf2-wait-%j.out
#SBATCH --error=mnf2-wait-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
root="$SLURM_SUBMIT_DIR/mnf2-wait-$SLURM_JOB_ID"
mkdir -p "$root"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
exe="$source_dir/build/src/xmvb-cpp.exe"
for variant in default passive; do
    mkdir -p "$root/$variant"
    cd "$root/$variant"
    cp "$source_dir/testdata/vbscf/MnF2.xmi" input.xmi
    sha256sum "$exe" input.xmi > hashes.txt
    launch=(env -u OMP_WAIT_POLICY)
    [[ "$variant" != passive ]] || launch=(env OMP_WAIT_POLICY=PASSIVE)
    srun --exact --cpus-per-task=32 --cpu-bind=cores \
        "${launch[@]}" XMVB_PHASE_TIMING=1 \
        /usr/bin/time -v -o resource.time "$exe" input.xmi \
        --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        > output.out 2> run.err
done
