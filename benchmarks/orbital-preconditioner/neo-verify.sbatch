#!/usr/bin/env bash
#SBATCH --job-name=neo-verify
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=15-00:00:00
#SBATCH --output=neo-verify-%j.out
#SBATCH --error=neo-verify-%j.err
set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export PATH="$dependencies/bin:$PATH" CONDA_PREFIX="$dependencies"
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1 LC_ALL=C
root="$SLURM_SUBMIT_DIR/neo-verify-$SLURM_JOB_ID"
mkdir -p "$root"
hostname > "$root/node.txt"
tests=(test_neo test_neo_response_solver test_spectral_trust_region
       test_active_two_electron_sparse test_ri_ao_h1e_hvp test_orthogonal_direct_ci)
cmake --build "$source_dir/build" --target xmvb "${tests[@]}" -j16 > "$root/build.log" 2>&1
for threads in 1 8 32; do
    for test in "${tests[@]}"; do
        srun --exact --cpus-per-task="$threads" --cpu-bind=cores \
            env OMP_NUM_THREADS="$threads" "$source_dir/build/src/$test" \
            > "$root/$test-$threads.out" 2>&1
    done
done
sha256sum "$source_dir/build/src/xmvb-cpp.exe" > "$root/hash.txt"
