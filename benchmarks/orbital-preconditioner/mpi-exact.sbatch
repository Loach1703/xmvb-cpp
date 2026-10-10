#!/usr/bin/env bash
#SBATCH --job-name=exact-mpi-paths
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15-00:00:00
#SBATCH --output=mpi-exact-%j.out
#SBATCH --error=mpi-exact-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06 openmpi/5.0.9
set -u
source_dir="$1"
baseline="$2"
previous="$3"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
root="$SLURM_SUBMIT_DIR/mpi-exact-$SLURM_JOB_ID"
mkdir -p "$root"
cp "$source_dir/build/src/xmvb-cpp.exe" "$root/current.exe"
cp "$previous" "$root/previous.exe"
cp "$baseline" "$root/baseline.exe"
anchor=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n 1)
cd "$root"
scontrol show job "$SLURM_JOB_ID" > allocation.txt
srun --exact --nodes=2 --ntasks=2 --ntasks-per-node=1 --cpus-per-task=1 \
    bash -c 'hostname > "node-${SLURM_NODEID}.txt"; lscpu > "cpu-${SLURM_NODEID}.txt"'
sha256sum "$root"/*.exe "$source_dir/testdata/vbscf/MnF2.xmi" \
    "$SLURM_SUBMIT_DIR/F2-auto.xmi" > hashes.txt
for molecule in MnF2 F2-auto; do
    for layout in baseline 1x32 2x16 2x32 previous-2x32; do
        [[ "$molecule" != F2-auto || "$layout" != previous-2x32 ]] || continue
        mkdir -p "$root/$molecule/$layout"
        cd "$root/$molecule/$layout"
        input="$source_dir/testdata/vbscf/MnF2.xmi"
        [[ "$molecule" != F2-auto ]] || input="$SLURM_SUBMIT_DIR/F2-auto.xmi"
        cp "$input" input.xmi
        ranks=1; threads=32; exe="$root/current.exe"
        [[ "$layout" != baseline ]] || exe="$root/baseline.exe"
        if [[ "$layout" == 2x* ]]; then ranks=2; threads="${layout#*x}"; fi
        if [[ "$layout" == previous-2x32 ]]; then ranks=2; exe="$root/previous.exe"; fi
        node_args=()
        [[ "$ranks" != 1 ]] || node_args=(--nodelist="$anchor")
        srun --exact --mpi=pmix_v5 --kill-on-bad-exit=1 --nodes="$ranks" --ntasks="$ranks" \
            --ntasks-per-node=1 --cpus-per-task="$threads" --cpu-bind=cores "${node_args[@]}" \
            env OMP_NUM_THREADS="$threads" XMVB_PHASE_TIMING=1 \
            bash -c 'exec /usr/bin/time -v -o "rank-${SLURM_PROCID}.time" "$1" input.xmi --optimizer-backend neo --gradient-tolerance 1e-3 --energy-tolerance 1e-7' \
            _ "$exe" > output.out 2> run.err
    done
done
