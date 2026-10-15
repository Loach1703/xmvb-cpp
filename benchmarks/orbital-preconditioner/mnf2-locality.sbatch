#!/usr/bin/env bash
#SBATCH --job-name=mnf2-locality
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15-00:00:00
#SBATCH --output=mnf2-locality-%j.out
#SBATCH --error=mnf2-locality-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export PATH="$dependencies/bin:$PATH" LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export LC_ALL=C OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
unset OMP_WAIT_POLICY GOMP_SPINCOUNT
root="$SLURM_SUBMIT_DIR/mnf2-locality-$SLURM_JOB_ID"
mkdir -p "$root"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
scontrol show job "$SLURM_JOB_ID" > "$root/allocation.txt"
cp "$source_dir/build/src/xmvb-cpp.exe" "$root/current.exe"
for variant in prefix-16 compact-16 compact-8 compact-32 compact-1; do
    threads="${variant##*-}"
    placement=(--compact)
    [[ "$variant" != prefix-16 ]] || placement=(--nested)
    mkdir -p "$root/$variant"
    cd "$root/$variant"
    cp "$source_dir/testdata/vbscf/MnF2.xmi" input.xmi
    sha256sum "$root/current.exe" input.xmi > hashes.txt
    srun --exact --cpus-per-task=32 --cpu-bind=verbose,cores \
        env OMP_NUM_THREADS="$threads" OMP_DISPLAY_ENV=VERBOSE OMP_DISPLAY_AFFINITY=TRUE XMVB_PHASE_TIMING=1 \
        python3 "$source_dir/benchmarks/orbital-preconditioner/affinity_probe.py" \
        --threads "$threads" "${placement[@]}" --record affinity.json -- \
        /usr/bin/time -v -o resource.time "$root/current.exe" input.xmi \
        --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        > output.out 2> run.err
done
