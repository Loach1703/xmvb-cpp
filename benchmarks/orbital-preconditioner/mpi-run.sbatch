#!/usr/bin/env bash
#SBATCH --job-name=mnf2-mpi-shards
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15-00:00:00
#SBATCH --output=mpi-run-%j.out
#SBATCH --error=mpi-run-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06 openmpi/5.0.9
set -u
source_dir="$1"
baseline="$2"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export PATH="$dependencies/bin:$PATH"
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
root="$SLURM_SUBMIT_DIR/mpi-run-$SLURM_JOB_ID"
mkdir -p "$root"
cd "$root"
cp "$source_dir/build/src/xmvb-cpp.exe" "$root/xmvb-cpp.exe"
cp "$source_dir/build/src/check_integral_shards" "$root/check_integral_shards.exe"
ompi_info --version > mpi-version.txt
scontrol show job "$SLURM_JOB_ID" > allocation.txt
anchor=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n 1)
srun --exact --nodes=2 --ntasks=2 --ntasks-per-node=1 --cpus-per-task=1 \
    bash -c 'hostname > "node-${SLURM_NODEID}.txt"; lscpu > "cpu-${SLURM_NODEID}.txt"'
sha256sum "$baseline" "$root/xmvb-cpp.exe" \
    "$root/check_integral_shards.exe" "$source_dir/testdata/vbscf/MnF2.xmi" > hashes.txt
if [[ -n ${3:-} ]]; then
    # Reuse the completed native pair when only diagnostic launch plumbing
    # changes; do not rerun expensive SCF work to obtain redundant timings.
    cp -a "$3/baseline" "$3/updated" "$root/"
    printf '%s\n' "$3" > "$root/native-reference.txt"
else
for variant in baseline updated; do
    mkdir -p "$root/$variant"
    cd "$root/$variant"
    cp "$source_dir/testdata/vbscf/MnF2.xmi" input.xmi
    exe="$baseline"
    [[ "$variant" != updated ]] || exe="$root/xmvb-cpp.exe"
    srun --exact --nodes=1 --nodelist="$anchor" --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
        env OMP_NUM_THREADS=32 XMVB_PHASE_TIMING=1 \
        /usr/bin/time -v -o resource.time "$exe" input.xmi \
        --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        > output.out 2> run.err
done
fi
cd "$root"
# Keep total cores equal for the first comparison, then measure 64-core scaling.
for layout in 1x32 2x16 2x32; do
    ranks="${layout%x*}"
    threads="${layout#*x}"
    mkdir -p "$root/$layout"
    cd "$root/$layout"
    node_args=()
    [[ "$ranks" != 1 ]] || node_args=(--nodelist="$anchor")
    srun --exact --mpi=pmix_v5 --nodes="$ranks" --ntasks="$ranks" \
        "${node_args[@]}" \
        --ntasks-per-node=1 --cpus-per-task="$threads" --cpu-bind=cores \
        env OMP_NUM_THREADS="$threads" \
        bash -c 'exec /usr/bin/time -v -o "rank-${SLURM_PROCID}.time" "$1" "$2"' \
        _ "$root/check_integral_shards.exe" "$source_dir/testdata/vbscf/MnF2.xmi" \
        > output.out 2> run.err
done
cd "$root"
python3 "$source_dir/benchmarks/orbital-preconditioner/mpi-check.py" "$root"
