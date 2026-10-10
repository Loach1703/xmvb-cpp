#!/usr/bin/env bash
#SBATCH --job-name=mnf2-mpi-scf
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15-00:00:00
#SBATCH --output=mpi-scf-%j.out
#SBATCH --error=mpi-scf-%j.err
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
root="$SLURM_SUBMIT_DIR/mpi-scf-$SLURM_JOB_ID"
mkdir -p "$root"
cp "$source_dir/build/src/xmvb-cpp.exe" "$root/xmvb-cpp.exe"
anchor=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n 1)
cd "$root"
scontrol show job "$SLURM_JOB_ID" > allocation.txt
srun --exact --nodes=2 --ntasks=2 --ntasks-per-node=1 --cpus-per-task=1 \
    bash -c 'hostname > "node-${SLURM_NODEID}.txt"; lscpu > "cpu-${SLURM_NODEID}.txt"'
sha256sum "$baseline" "$root/xmvb-cpp.exe" "$source_dir/testdata/vbscf/MnF2.xmi" > hashes.txt
for layout in baseline 1x32 2x16 2x32; do
    mkdir -p "$root/$layout"
    cd "$root/$layout"
    cp "$source_dir/testdata/vbscf/MnF2.xmi" input.xmi
    ranks=1
    threads=32
    exe="$root/xmvb-cpp.exe"
    [[ "$layout" != baseline ]] || exe="$baseline"
    if [[ "$layout" == 2x* ]]; then ranks=2; threads="${layout#*x}"; fi
    node_args=()
    [[ "$ranks" != 1 ]] || node_args=(--nodelist="$anchor")
    srun --exact --mpi=pmix_v5 --kill-on-bad-exit=1 --nodes="$ranks" --ntasks="$ranks" \
        --ntasks-per-node=1 --cpus-per-task="$threads" --cpu-bind=cores "${node_args[@]}" \
        env OMP_NUM_THREADS="$threads" XMVB_PHASE_TIMING=1 \
        bash -c 'exec /usr/bin/time -v -o "rank-${SLURM_PROCID}.time" "$1" input.xmi --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7' \
        _ "$exe" > output.out 2> run.err
done
srun --exact --mpi=pmix_v5 --kill-on-bad-exit=1 --nodes=2 --ntasks=2 \
    --ntasks-per-node=1 --cpus-per-task=16 --cpu-bind=cores \
    env OMP_NUM_THREADS=16 "$source_dir/build/src/check_integral_service" \
    "$source_dir/testdata/vbscf/F2.xmi" > "$root/service-check.out" 2> "$root/service-check.err"
