#!/usr/bin/env bash
#SBATCH --job-name=fe-porphine-optimizers
#SBATCH --partition=CPU-256C768GB
#SBATCH --qos=qos_cpu-256c768gb
#SBATCH --account=taoxia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -eo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || exit 64
. /etc/profile
module purge
module load anaconda3/25.06
set -u
source_dir="$1"
dependencies=/home/guqqgroup/taoxia/software/conda-envs/xmvb-cpp-89c74a6
export LC_ALL=C OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 GOTO_NUM_THREADS=1
export CCACHE_DIR=/home/guqqgroup/taoxia/.cache/ccache-xmvb-cpp
hostname > node.txt
lscpu > cpu.txt
module list > modules.txt 2>&1
env | sort | rg '^(OMP_|OPENBLAS_|MKL_|BLIS_|GOTO_|SLURM_)' > allocation.txt
(
    export PATH="$dependencies/bin:$PATH" CONDA_PREFIX="$dependencies"
    export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
    cmake -S "$source_dir" -B "$source_dir/build" -G Ninja \
        -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
        -DCMAKE_PREFIX_PATH="$dependencies" \
        -DCMAKE_C_COMPILER="$dependencies/bin/cc" \
        -DCMAKE_CXX_COMPILER="$dependencies/bin/c++" \
        -DCMAKE_Fortran_COMPILER="$dependencies/bin/gfortran" \
        -DXMVB_ENABLE_EIGEN_BLAS_BACKEND=ON
    cmake --build "$source_dir/build" --target xmvb -j16
) > build.log 2>&1
srun --cpu-bind=cores /usr/bin/time -v -o prepare.time \
    /home/guqqgroup/taoxia/conda-envs/pyscf/bin/python -u prepare.py "$source_dir" \
    > prepare.out 2> prepare.err
export PATH="$dependencies/bin:$PATH" CONDA_PREFIX="$dependencies"
export LD_LIBRARY_PATH="$dependencies/lib:${LD_LIBRARY_PATH:-}"
exe="$source_dir/build/src/xmvb-cpp.exe"
sha256sum "$exe" > binary.sha256
printf 'state,iscf,exit_code\n' > exit_codes.csv
for state in singlet triplet quintet; do
    for iscf in 5 8 9; do
        label="$state-iscf$iscf"
        mkdir -p "$label"
        cp "$label.xmi" "$label/input.xmi"
        sha256sum "$label.xmi" >> inputs.sha256
        printf '%s START %s\n' "$(date -Is)" "$label" >> progress.txt
        status=0
        (
            cd "$label"
            # The safety bound is wall time, not a 30-step truncation. A timed-out
            # calculation remains an incomplete result, never a speedup claim.
            srun --cpu-bind=cores /usr/bin/time -v -o resource.time \
                timeout --signal=TERM --kill-after=60s 2h "$exe" input.xmi \
                --gradient-tolerance 1e-3 --energy-tolerance 1e-7 \
                > output.out 2> run.err
        ) || status=$?
        printf '%s,%s,%s\n' "$state" "$iscf" "$status" >> exit_codes.csv
        printf '%s END %s %s\n' "$(date -Is)" "$label" "$status" >> progress.txt
    done
done
date -Is > finished.txt
