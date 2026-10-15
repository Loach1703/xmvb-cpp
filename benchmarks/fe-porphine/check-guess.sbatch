#!/usr/bin/env bash
#SBATCH --job-name=fe-guess-check
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=01:00:00
#SBATCH --output=check-guess-%j.out
#SBATCH --error=check-guess-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
checkpoint="$2"
binary="$3"
run_dir="$SLURM_SUBMIT_DIR/guess-check-$SLURM_JOB_ID"
mkdir -p "$run_dir"
cp "$SLURM_SUBMIT_DIR/geometry.xyz" "$SLURM_SUBMIT_DIR/prepare.py" \
    "$SLURM_SUBMIT_DIR/audit_guess.py" "$SLURM_SUBMIT_DIR/check_guess.py" "$run_dir/"
cd "$run_dir"
hostname > node.txt
sha256sum "$checkpoint" "$binary" prepare.py check_guess.py > hashes.txt
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_MAX_ACTIVE_LEVELS=1
python_bin=/home/guqqgroup/taoxia/conda-envs/pyscf/bin/python
srun --cpu-bind=cores /usr/bin/time -v -o prepare.time \
    "$python_bin" -u prepare.py "$source_dir" "$checkpoint" > prepare.out 2> prepare.err
srun --cpu-bind=cores "$python_bin" audit_guess.py "$run_dir" > audit.json
srun --cpu-bind=cores /usr/bin/time -v -o casci.time \
    "$python_bin" -u check_guess.py > casci.out 2> casci.err
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
# One step only: record the initial CI energy and one update, not convergence.
srun --cpu-bind=cores /usr/bin/time -v -o cpp.time \
    "$binary" quintet-iscf8.xmi --max-iterations 1 > cpp.out 2> cpp.err
