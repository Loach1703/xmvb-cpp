#!/usr/bin/env bash
#SBATCH --job-name=241-neo-threads
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --mem=96G
#SBATCH --time=15-00:00:00
#SBATCH --output=thread-scaling-%j.out
#SBATCH --error=thread-scaling-%j.err
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
exe="$source_dir/build/src/xmvb-cpp.exe"
vtune=/opt/oneapi/2025.2.1.44/vtune/2025.4/bin64/vtune
[[ -x "$exe" && -x "$vtune" ]] || exit 66
root="$SLURM_SUBMIT_DIR/thread-scaling-$SLURM_JOB_ID"
mkdir -p "$root"
hostname > "$root/node.txt"
lscpu > "$root/cpu.txt"
module list > "$root/modules.txt" 2>&1
env | sort > "$root/environment.txt"
cp "$SLURM_SUBMIT_DIR/thread-scaling.sbatch" "$root/"
for threads in 1 8 32; do
    run_dir="$root/native-$threads"
    mkdir -p "$run_dir"
    cd "$run_dir"
    cp "$source_dir/testdata/vbscf/241_vb.xmi" input.xmi
    sha256sum "$exe" input.xmi > hashes.txt
    export OMP_NUM_THREADS="$threads"
    srun --exact --ntasks=1 --cpus-per-task="$threads" --cpu-bind=cores \
        env XMVB_PHASE_TIMING=1 /usr/bin/time -v -o resource.time \
        "$exe" input.xmi --optimizer-backend neo \
        --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
        > output.out 2> run.err
done
mkdir -p "$root/vtune"
cd "$root/vtune"
cp "$source_dir/testdata/vbscf/241_vb.xmi" input.xmi
export OMP_NUM_THREADS=32
sha256sum "$exe" input.xmi > hashes.txt
srun --exact --ntasks=1 --cpus-per-task=32 --cpu-bind=cores \
    "$vtune" -collect hotspots -knob sampling-mode=sw -result-dir capture -- \
    "$exe" input.xmi --optimizer-backend neo \
    --gradient-tolerance 1e-3 --energy-tolerance 1e-7 > output.out 2> run.err
srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
    "$vtune" -report summary -r capture > summary.txt 2> summary.err
srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
    "$vtune" -report hotspots -r capture -group-by function -format csv \
    > hotspots.csv 2> hotspots.err
srun --exact --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
    "$vtune" -report hotspots -r capture -group-by thread -format csv \
    > threads.csv 2> threads.err
